"""N-glycosylation feature extraction + prediction (StackGlyEmbed), 100% local/offline.

Adaptation of ``extractFeatures.py`` + ``predict.py`` from the original repo
(github.com/nafcoder/StackGlyEmbed) for use as a subprocess-invoked
script from ``ProtStackGlyEmbedPrediction``. This script is vendored inside
this plugin's own package (not the user's separately-cloned StackGlyEmbed
repo) because the classifier pickles (``power_transformer_*.sav``,
``base_layer_pickle_files/``) live inside that clone and are referenced
only via ``--models-dir`` (never hardcoded), while this extraction/
prediction logic is original integration code, kept in this plugin's own
versioned tree.

Changes vs. the original:

1. **ESM-2 via ``transformers.EsmModel`` instead of ``torch.hub.load``**: the
   original calls ``torch.hub.load("facebookresearch/esm:main", ...)``, which
   hits the network on every run. ``facebook/esm2_t33_650M_UR50D`` via
   ``transformers`` is the same model (same author, same checkpoint),
   supports fully offline loading once cached (``HF_HUB_OFFLINE=1``, see
   below), and produces the same per-layer representations (33 layers, dim
   1280) the original code consumes.

2. **ProtT5 points at a LOCAL path** (by default, the weights already
   downloaded for TMbed, ``Rostlab/prot_t5_xl_half_uniref50-enc`` -- same
   encoder as the original's ``Rostlab/prot_t5_xl_uniref50``, just fp16/no
   decoder) instead of the remote ID the original always resolved against
   the Hub.

3. **ProteinBERT points at a plugin-local dir** (``--proteinbert-dir``,
   defaults to the library's own ``~/proteinbert_models`` only if unset) --
   ``download_model_dump_if_not_exists=False`` is passed so it *fails*
   with a clear error instead of silently attempting network access
   mid-run if the dump is missing for any reason (it is downloaded once,
   during plugin installation, not here).

4. **Input/output/model paths via CLI** instead of cwd-relative hardcoded
   filenames (``dataset.txt``, ``features.csv``, ``predicted_values.txt``,
   pickles in ``./base_layer_pickle_files/`` in the original): allows
   multiple concurrent/successive invocations without clobbering each other
   or depending on a specific external clone's directory as cwd.
   ``--dataset`` and ``--output-dir`` are confined to the working directory
   (see :func:`resolveInsideWorkingDir`), which under Scipion is the project
   directory the protocol's ``runJob`` passes as cwd.

5. **ProteinBERT's global embedding is computed ONCE per protein** (outside
   the site loop), not once per candidate site as in the original
   (``getRepresentation`` was called inside the ``site_position`` loop,
   recomputing the same deterministic forward pass N times for N sites of
   the same protein). Numerically identical result, just avoids redundant
   work.

The ``dataset.txt`` format (``Protein_id,site_1,site_2,...`` followed by the
sequence on the next line) and the column order/dimensionality of
``features.csv`` (ProteinBERT global + ESM-2 windowed mean + ProtT5
residue-point, in that order) are kept EXACT relative to the original: the
classifiers in ``base_layer_pickle_files/`` and the ``power_transformer_*.sav``
files were trained against that precise order and are not touched here.
"""

import argparse
import os

os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

import pickle  # noqa: E402
import re  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402
from sklearn.preprocessing import PowerTransformer  # noqa: E402,F401 (needed to unpickle power_transformer_*.sav)
from tensorflow import keras  # noqa: E402
from transformers import AutoTokenizer, EsmModel, T5EncoderModel, T5Tokenizer  # noqa: E402

from proteinbert import load_pretrained_model  # noqa: E402

_WINDOW_SIZE = 15
_DEFAULT_ESM_MODEL = "facebook/esm2_t33_650M_UR50D"


def getModelWithGlobalEmbeddingAsOutputs(model):
    """Rebuilds the ProteinBERT model to expose the global embedding (see original README/extractFeatures.py)."""
    globalLayers = [
        layer.output
        for layer in model.layers
        if len(layer.output.shape) == 2 and layer.name in ["global-merge2-norm-block6"]
    ]
    concatenated = keras.layers.Concatenate(name="last-Window-layers")(globalLayers)
    return keras.models.Model(inputs=model.inputs, outputs=concatenated)


def getProteinbertRepresentation(pretrainedModelGenerator, inputEncoder, seq: str) -> np.ndarray:
    encodedX = inputEncoder.encode_X([seq], len(seq) + 2)
    model = getModelWithGlobalEmbeddingAsOutputs(pretrainedModelGenerator.create_model(len(seq) + 2))
    return np.array(model.predict(encodedX, batch_size=2))[0]


def getEsm2Embedding(tokenizer, model, seq: str) -> np.ndarray:
    """Per-residue ESM-2 embedding (last-layer representations, no CLS/EOS)."""
    chunks = [seq[i : i + 1024] for i in range(0, len(seq), 1024)]
    final = np.zeros((1, model.config.hidden_size))
    for chunk in chunks:
        tokens = tokenizer(chunk, return_tensors="pt")
        with torch.no_grad():
            out = model(**tokens)
        rep = out.last_hidden_state[0, 1:-1].numpy()
        final = np.concatenate((final, rep), axis=0)
    return np.delete(final, 0, axis=0)


def getProtT5Embedding(tokenizer, model, seq: str) -> np.ndarray:
    """Per-residue ProtT5 embedding (same 8797 aa chunking as the original script)."""
    chunks = [seq[i : i + 8797] for i in range(0, len(seq), 8797)]
    final = np.zeros((1, model.config.d_model))
    for chunk in chunks:
        spaced = " ".join(list(re.sub(r"[UZOB]", "X", chunk)))
        ids = tokenizer([spaced], add_special_tokens=True, padding="longest", return_tensors="pt")
        with torch.no_grad():
            out = model(input_ids=ids["input_ids"], attention_mask=ids["attention_mask"])
        emb = out.last_hidden_state[0, : len(chunk)].numpy()
        final = np.concatenate((final, emb), axis=0)
    return np.delete(final, 0, axis=0)


def resolveInsideWorkingDir(rawPath: Path, argName: str) -> Path:
    """Canonicalizes a path given on the command line and keeps it inside the cwd.

    ``--dataset`` and ``--output-dir`` are the only two paths this script reads
    from / writes to on its caller's behalf. Under Scipion both arrive relative
    to the project directory, which is the cwd ``runJob`` hands over (e.g.
    ``Runs/<runId>_ProtStackGlyEmbedPrediction/extra``), so anchoring them there
    matches the real caller exactly while stopping a ``..`` sequence from
    reaching ``mkdir()`` or ``read_text()`` outside the project.

    The weight/model paths (``--models-dir``, ``--proteinbert-dir``,
    ``--t5-model-path``) are plugin installation paths, not per-run caller
    input, and live outside the project on purpose: they are NOT confined here.
    """
    root = os.path.realpath(os.getcwd())
    resolved = os.path.realpath(os.path.join(root, str(rawPath)))
    if resolved != root and not resolved.startswith(root + os.sep):
        raise ValueError(
            f"{argName} must stay inside the working directory '{root}', got '{rawPath}'"
        )
    return Path(resolved)


def extractFeatures(datasetPath: Path, outputDir: Path, t5ModelPath: str, esmModelName: str,
                    proteinbertDir: str = None) -> Path:
    """Generates ``features.csv`` (ProteinBERT + ESM-2 + ProtT5) for each site in ``dataset.txt``."""
    outputDir.mkdir(parents=True, exist_ok=True)

    proteinbertKwargs = {"download_model_dump_if_not_exists": False}
    if proteinbertDir:
        proteinbertKwargs["local_model_dump_dir"] = proteinbertDir
    print(f"Loading ProteinBERT (local, {proteinbertDir or '~/proteinbert_models'}/default.pkl)...", flush=True)
    pretrainedModelGenerator, inputEncoder = load_pretrained_model(**proteinbertKwargs)

    print(f"Loading ESM-2 650M ({esmModelName}, offline local)...", flush=True)
    esmTokenizer = AutoTokenizer.from_pretrained(esmModelName)
    esmModel = EsmModel.from_pretrained(esmModelName).eval()

    print(f"Loading ProtT5 ({t5ModelPath}, offline local)...", flush=True)
    t5Tokenizer = T5Tokenizer.from_pretrained(t5ModelPath, do_lower_case=False)
    t5Model = T5EncoderModel.from_pretrained(t5ModelPath).eval()

    lines = datasetPath.read_text().splitlines()

    proteinbertRows, esmRows, t5Rows = [], [], []
    for i in range(0, len(lines), 2):
        header = lines[i].split(",")
        proteinId = header[0]
        positions = [int(p) for p in header[1:]]
        seq = lines[i + 1]

        print(f"[{proteinId}] {len(positions)} candidate site(s), {len(seq)} aa", flush=True)
        pbFull = getProteinbertRepresentation(pretrainedModelGenerator, inputEncoder, seq)
        esmFull = getEsm2Embedding(esmTokenizer, esmModel, seq)
        t5Full = getProtT5Embedding(t5Tokenizer, t5Model, seq)

        for pos in positions:
            proteinbertRows.append(pbFull)
            start = max(pos - _WINDOW_SIZE - 1, 0)
            end = min(pos + _WINDOW_SIZE, len(seq))
            esmRows.append(np.mean(esmFull[start:end, :], axis=0))
            t5Rows.append(t5Full[pos - 1])

    features = np.concatenate([np.array(proteinbertRows), np.array(esmRows), np.array(t5Rows)], axis=1)
    featuresPath = outputDir / "features.csv"
    np.savetxt(featuresPath, features, delimiter=",")
    return featuresPath


def preprocess(featureX: np.ndarray, stage: int, modelsDir: Path) -> np.ndarray:
    with open(modelsDir / f"power_transformer_{stage}.sav", "rb") as f:
        pt = pickle.load(f)
    return pt.transform(featureX)


def baseLayerPredictions(featureX: np.ndarray, modelsDir: Path) -> np.ndarray:
    testX = preprocess(featureX, 2, modelsDir)
    total = np.zeros((len(testX), 1), dtype=float)
    pickleDir = modelsDir / "base_layer_pickle_files"

    for i in range(10):
        for baseClassifier in ("SVM", "XGB", "KNN"):
            with open(pickleDir / f"{baseClassifier}_base_layer_{i}.sav", "rb") as f:
                model = pickle.load(f)
            yProba = model.predict_proba(testX)[:, 1].reshape(-1, 1)
            total = np.concatenate((total, yProba), axis=1)

    return np.delete(total, 0, axis=1)


def predict(featuresPath: Path, outputDir: Path, modelsDir: Path) -> Path:
    """Applies the already-trained classifier stack (base layer + meta-SVM), unmodified."""
    featureX = np.loadtxt(featuresPath, delimiter=",")
    if featureX.ndim == 1:
        featureX = featureX.reshape(1, -1)

    x = preprocess(featureX, 1, modelsDir)
    blp = baseLayerPredictions(x, modelsDir)
    x = np.concatenate((x, blp), axis=1)
    x = preprocess(x, 3, modelsDir)

    with open(modelsDir / "base_layer_pickle_files" / "SVM_meta_layer.sav", "rb") as f:
        clf = pickle.load(f)
    yPred = clf.predict(x)
    yProba = clf.predict_proba(x)[:, 1]

    predictedPath = outputDir / "predicted_values.csv"
    np.savetxt(predictedPath, np.column_stack([yPred, yProba]), delimiter=",", fmt="%.6f",
               header="prediction,probability", comments="")
    return predictedPath


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path, help="dataset.txt (original StackGlyEmbed format)")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to write features.csv and predicted_values.csv")
    parser.add_argument("--models-dir", required=True, type=Path,
                         help="'prediction/' folder of the StackGlyEmbed clone (power_transformer_*.sav + base_layer_pickle_files/)")
    parser.add_argument("--t5-model-path", required=True, help="Local path to ProtT5 weights")
    parser.add_argument("--esm-model-name", default=_DEFAULT_ESM_MODEL, help="HF Hub ID of the ESM-2 model (offline if already cached)")
    parser.add_argument("--proteinbert-dir", default=None,
                         help="Local dir containing ProteinBERT's default.pkl dump (defaults to ~/proteinbert_models if unset)")
    args = parser.parse_args()

    datasetPath = resolveInsideWorkingDir(args.dataset, "--dataset")
    outputDir = resolveInsideWorkingDir(args.output_dir, "--output-dir")

    featuresPath = extractFeatures(datasetPath, outputDir, args.t5_model_path, args.esm_model_name,
                                   args.proteinbert_dir)
    predictedPath = predict(featuresPath, outputDir, args.models_dir)
    print(f"-> Predictions saved to: {predictedPath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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


def _get_model_with_global_embedding_as_outputs(model):
    """Rebuilds the ProteinBERT model to expose the global embedding (see original README/extractFeatures.py)."""
    global_layers = [
        layer.output
        for layer in model.layers
        if len(layer.output.shape) == 2 and layer.name in ["global-merge2-norm-block6"]
    ]
    concatenated = keras.layers.Concatenate(name="last-Window-layers")(global_layers)
    return keras.models.Model(inputs=model.inputs, outputs=concatenated)


def _get_proteinbert_representation(pretrained_model_generator, input_encoder, seq: str) -> np.ndarray:
    encoded_x = input_encoder.encode_X([seq], len(seq) + 2)
    model = _get_model_with_global_embedding_as_outputs(pretrained_model_generator.create_model(len(seq) + 2))
    return np.array(model.predict(encoded_x, batch_size=2))[0]


def _get_esm2_embedding(tokenizer, model, seq: str) -> np.ndarray:
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


def _get_prott5_embedding(tokenizer, model, seq: str) -> np.ndarray:
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


def extract_features(dataset_path: Path, output_dir: Path, t5_model_path: str, esm_model_name: str,
                      proteinbert_dir: str = None) -> Path:
    """Generates ``features.csv`` (ProteinBERT + ESM-2 + ProtT5) for each site in ``dataset.txt``."""
    output_dir.mkdir(parents=True, exist_ok=True)

    proteinbertKwargs = {"download_model_dump_if_not_exists": False}
    if proteinbert_dir:
        proteinbertKwargs["local_model_dump_dir"] = proteinbert_dir
    print(f"Loading ProteinBERT (local, {proteinbert_dir or '~/proteinbert_models'}/default.pkl)...", flush=True)
    pretrained_model_generator, input_encoder = load_pretrained_model(**proteinbertKwargs)

    print(f"Loading ESM-2 650M ({esm_model_name}, offline local)...", flush=True)
    esm_tokenizer = AutoTokenizer.from_pretrained(esm_model_name)
    esm_model = EsmModel.from_pretrained(esm_model_name).eval()

    print(f"Loading ProtT5 ({t5_model_path}, offline local)...", flush=True)
    t5_tokenizer = T5Tokenizer.from_pretrained(t5_model_path, do_lower_case=False)
    t5_model = T5EncoderModel.from_pretrained(t5_model_path).eval()

    lines = dataset_path.read_text().splitlines()

    proteinbert_rows, esm_rows, t5_rows = [], [], []
    for i in range(0, len(lines), 2):
        header = lines[i].split(",")
        protein_id = header[0]
        positions = [int(p) for p in header[1:]]
        seq = lines[i + 1]

        print(f"[{protein_id}] {len(positions)} candidate site(s), {len(seq)} aa", flush=True)
        pb_full = _get_proteinbert_representation(pretrained_model_generator, input_encoder, seq)
        esm_full = _get_esm2_embedding(esm_tokenizer, esm_model, seq)
        t5_full = _get_prott5_embedding(t5_tokenizer, t5_model, seq)

        for pos in positions:
            proteinbert_rows.append(pb_full)
            start = max(pos - _WINDOW_SIZE - 1, 0)
            end = min(pos + _WINDOW_SIZE, len(seq))
            esm_rows.append(np.mean(esm_full[start:end, :], axis=0))
            t5_rows.append(t5_full[pos - 1])

    features = np.concatenate([np.array(proteinbert_rows), np.array(esm_rows), np.array(t5_rows)], axis=1)
    features_path = output_dir / "features.csv"
    np.savetxt(features_path, features, delimiter=",")
    return features_path


def _preprocess(feature_x: np.ndarray, stage: int, models_dir: Path) -> np.ndarray:
    with open(models_dir / f"power_transformer_{stage}.sav", "rb") as f:
        pt = pickle.load(f)
    return pt.transform(feature_x)


def _base_layer_predictions(feature_x: np.ndarray, models_dir: Path) -> np.ndarray:
    test_x = _preprocess(feature_x, 2, models_dir)
    total = np.zeros((len(test_x), 1), dtype=float)
    pickle_dir = models_dir / "base_layer_pickle_files"

    for i in range(10):
        for base_classifier in ("SVM", "XGB", "KNN"):
            with open(pickle_dir / f"{base_classifier}_base_layer_{i}.sav", "rb") as f:
                model = pickle.load(f)
            y_proba = model.predict_proba(test_x)[:, 1].reshape(-1, 1)
            total = np.concatenate((total, y_proba), axis=1)

    return np.delete(total, 0, axis=1)


def predict(features_path: Path, output_dir: Path, models_dir: Path) -> Path:
    """Applies the already-trained classifier stack (base layer + meta-SVM), unmodified."""
    feature_x = np.loadtxt(features_path, delimiter=",")
    if feature_x.ndim == 1:
        feature_x = feature_x.reshape(1, -1)

    x = _preprocess(feature_x, 1, models_dir)
    blp = _base_layer_predictions(x, models_dir)
    x = np.concatenate((x, blp), axis=1)
    x = _preprocess(x, 3, models_dir)

    with open(models_dir / "base_layer_pickle_files" / "SVM_meta_layer.sav", "rb") as f:
        clf = pickle.load(f)
    y_pred = clf.predict(x)
    y_proba = clf.predict_proba(x)[:, 1]

    predicted_path = output_dir / "predicted_values.csv"
    np.savetxt(predicted_path, np.column_stack([y_pred, y_proba]), delimiter=",", fmt="%.6f",
               header="prediction,probability", comments="")
    return predicted_path


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

    features_path = extract_features(args.dataset, args.output_dir, args.t5_model_path, args.esm_model_name,
                                      args.proteinbert_dir)
    predicted_path = predict(features_path, args.output_dir, args.models_dir)
    print(f"-> Predictions saved to: {predicted_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

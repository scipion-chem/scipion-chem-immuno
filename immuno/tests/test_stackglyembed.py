import re

from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from ..protocols import ProtStackGlyEmbedPrediction


class TestStackGlyEmbedPrediction(BaseTest):
    NAME = 'GP120_P03377'
    DESCRIPTION = 'ENV_HV1BR (GP120), UniProt P03377'
    AMINOACIDSSEQ = (
        'MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVW'
        'ATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLGNAT'
        'NTNSSNTNSSSGEMMMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYTLTSCNTSVIT'
        'QACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVV'
        'IRSANFTDNAKTIIVQLNQSVEINCTRPNNNTRKSIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNA'
        'TLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTE'
        'GSDTITLPCRIKQFINMWQEVGKAMYAPPISGQIRCSSNITGLLLTRDGGNNNNGSEIFRPGGGDMRDN'
        'WRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGARSMTLTVQARQLLSG'
        'IVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSN'
        'KSLEQIWNNMTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKI'
        'FIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSL'
        'ALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAV'
        'AEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL'
    )

    # 3 windows, each isolating exactly ONE real N-X-[S/T] sequon of this
    # sequence (confirmed via a real regex scan: 32 total sequons found,
    # these 3 are far enough from their neighbours to stay isolated inside
    # a +-10 aa window). Expected verdicts/scores come from a real local run
    # of the vendored predict_local.py script (not estimated): 2
    # 'Glycosylated' + 1 'Not glycosylated', a genuine mixed result.
    WINDOWS = [(80, 100), (745, 765), (811, 831)]
    # window start -> (isolated absolute sequon site, verdict, score)
    WINDOW_SITES = {80: 88, 745: 755, 811: 821}
    EXPECTED = {
        88: ('Glycosylated', 0.579834),
        755: ('Not glycosylated', 0.354970),
        821: ('Glycosylated', 0.657464),
    }

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

        cls.protSeedROIs = cls._runDefSeqROIs(cls.protImportSeq)
        cls._waitOutput(cls.protSeedROIs, 'outputROIs', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.DESCRIPTION,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    @classmethod
    def _runDefSeqROIs(cls, inProt):
        inROIs = '\n'.join(
            '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
                i, start, end, cls.AMINOACIDSSEQ[start - 1:end]
            )
            for i, (start, end) in enumerate(cls.WINDOWS, 1)
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(inProt)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')

        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def runStackGlyEmbed(self):
        protStackGly = self.newProtocol(ProtStackGlyEmbedPrediction)
        protStackGly.inputROIs.set(self.protSeedROIs)
        protStackGly.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protStackGly, wait=True)
        return protStackGly

    def test(self):
        protStackGly = self.runStackGlyEmbed()

        outROIs = getattr(protStackGly, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), len(self.WINDOWS))

        # NOTE: iterating a Scipion SetOfXXX reuses the same underlying
        # Python object per row (the sqlite cursor) -- .clone() each item
        # before storing it, or every value below silently ends up pointing
        # at the LAST row's state (the dict keys are still correct, since
        # they're evaluated synchronously during iteration; only the stored
        # object references are the problem).
        roisByStart = {roi.getROIIdx(): roi.clone() for roi in outROIs}
        for start, end in self.WINDOWS:
            site = self.WINDOW_SITES[start]
            expectedVerdict, expectedScore = self.EXPECTED[site]

            roi = roisByStart[start]
            summary = roi._glycoSequonSummary.get()
            match = re.search(rf'{site}:([^(]+)\(([\d.]+)\)', summary)
            self.assertIsNotNone(match, f'no entry for site {site} in summary {summary!r}')
            self.assertEqual(match.group(1), expectedVerdict)
            self.assertAlmostEqual(float(match.group(2)), expectedScore, places=3)
            self.assertEqual(roi._hasGlycoSequon.get(), expectedVerdict == 'Glycosylated')

import re

from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from ..protocols import ProtSignalPPrediction


class TestSignalPPrediction(BaseTest):
    NAME = 'SIGNALP_TEST_SEQ'
    DESCRIPTION = 'GP120 N-term fragment (no real signal peptide) + human insulin signal peptide (real SP)'
    PEPTIDES = [
        'MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEA',
        'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT',
    ]
    SPACER = 'GGG'
    AMINOACIDSSEQ = SPACER.join(PEPTIDES)

    # Real SignalP-6.0 output (organism=other, mode=slow-sequential), from
    # a direct local run of the real binary -- not estimated. Peptide 1
    # (GP120 N-term fragment) is a genuine negative control (no real
    # signal peptide); peptide 2 is the REAL human insulin signal peptide
    # (a well-known textbook example), confirming the positive-detection
    # path with a real predicted cleavage site. The CS position string's
    # embedded probability varies by +-0.0001 between real runs (minor
    # inference non-determinism confirmed empirically), so the test below
    # checks the "CS pos: X-Y" prefix exactly and the embedded probability
    # with a tolerance, rather than exact string equality.
    EXPECTED = {
        'MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEA': ('OTHER', 0.979668, 0.012823, None),
        'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT':
            ('SP', 0.000200, 0.999143, (24, 25, 0.9753)),
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
    def _getWindows(cls):
        windows = []
        cursor = 0
        for pep in cls.PEPTIDES:
            start = cls.AMINOACIDSSEQ.index(pep, cursor) + 1
            end = start + len(pep) - 1
            windows.append((start, end))
            cursor = end
        return windows

    @classmethod
    def _runDefSeqROIs(cls, inProt):
        windows = cls._getWindows()
        inROIs = '\n'.join(
            '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
                i, start, end, cls.AMINOACIDSSEQ[start - 1:end]
            )
            for i, (start, end) in enumerate(windows, 1)
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(inProt)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')

        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def test(self):
        protSignalP = self.newProtocol(ProtSignalPPrediction)
        protSignalP.inputROIs.set(self.protSeedROIs)
        protSignalP.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protSignalP, wait=True)

        outROIs = getattr(protSignalP, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), len(self.PEPTIDES))

        for roi in outROIs:
            seq = roi.getROISequence()
            expectedPrediction, expectedProbOther, expectedProbSp, expectedCs = self.EXPECTED[seq]
            self.assertEqual(roi._signalpPrediction.get(), expectedPrediction)
            self.assertAlmostEqual(roi._signalpProbOther.get(), expectedProbOther, places=4)
            self.assertAlmostEqual(roi._signalpProbSp.get(), expectedProbSp, places=4)

            csPosition = roi._signalpCsPosition.get()
            if expectedCs is None:
                self.assertEqual(csPosition, '')
            else:
                expectedStart, expectedEnd, expectedProb = expectedCs
                match = re.match(r'CS pos: (\d+)-(\d+)\. Pr: ([\d.]+)', csPosition)
                self.assertIsNotNone(match, f'unexpected CS position format: {csPosition!r}')
                self.assertEqual(int(match.group(1)), expectedStart)
                self.assertEqual(int(match.group(2)), expectedEnd)
                self.assertAlmostEqual(float(match.group(3)), expectedProb, places=2)

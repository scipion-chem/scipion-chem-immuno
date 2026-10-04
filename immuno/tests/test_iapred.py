from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from ..protocols import ProtIApredPrediction


class TestIApredPrediction(BaseTest):
    NAME = 'IAPRED_TEST_SEQ'
    DESCRIPTION = 'GP120 N-term fragment (48aa) + a real Flu M1 epitope (9aa, real <20aa case)'
    PEPTIDES = ['MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEA', 'GILGFVFTL']
    SPACER = 'GGG'
    AMINOACIDSSEQ = SPACER.join(PEPTIDES)

    # Real IApred output from a direct local run of the real binary -- not
    # estimated. The 9aa peptide genuinely triggers IApred's real <20aa
    # "Sequence too short" text-not-number limitation.
    EXPECTED = {
        'MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEA': ('Low', -0.92),
        'GILGFVFTL': ('Not evaluated (sequence < 20 aa)', None),
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
        protIApred = self.newProtocol(ProtIApredPrediction)
        protIApred.inputROIs.set(self.protSeedROIs)
        protIApred.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protIApred, wait=True)

        outROIs = getattr(protIApred, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), len(self.PEPTIDES))

        for roi in outROIs:
            seq = roi.getROISequence()
            expectedCategory, expectedScore = self.EXPECTED[seq]
            self.assertEqual(roi._iapredCategory.get(), expectedCategory)
            if expectedScore is None:
                self.assertIsNone(roi._iapredScore.get())
            else:
                self.assertAlmostEqual(roi._iapredScore.get(), expectedScore, places=2)

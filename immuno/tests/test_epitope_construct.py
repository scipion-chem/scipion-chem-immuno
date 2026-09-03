import os

from pwchem.objects import SetOfSequenceROIs
from pwchem.protocols import ProtDefineSeqROI
from pwem.protocols import ProtImportSequence
from pyworkflow.object import Boolean, Float, Integer, String
from pyworkflow.tests import BaseTest, setupTestProject

from ..protocols import ProtEpitopeConstructAssembly


class TestEpitopeConstructAssembly(BaseTest):
    # Real peptide substrings reused across this plugin's test fixtures
    # (HIV-1 Env GP120), not synthetic filler:
    BCELL_GOOD = 'QETLREHYQYVGKLAGRLKEASEGS'   # Non-Allergen, meanScore=0.8 -> ranked first
    BCELL_ALLERGEN = 'REHYQYVGKLAGRLKEASEG'    # Allergen, meanScore=0.2 -> not excluded
    # (see assembly.select_bcell_candidates docstring) -- ranked second, included
    HTL_A = 'WKNDMVEQM'   # nProm=5, minRank=0.5 -> ranked first
    HTL_B = 'IRIQRGPGR'   # nProm=3, minRank=1.0 -> ranked second
    CTL_A = 'RAIEAQQHL'   # netcleaveMatch=True -> ranked first despite worse nProm/minRank
    CTL_B = 'NAKTIIVQL'   # netcleaveMatch=False, nProm=6, minRank=0.3

    # Expected construct: B-cell (2, intra-linker KK, BCELL_GOOD before
    # BCELL_ALLERGEN by meanScore -- allergenicity no longer excludes, see
    # above) -> inter-block GPGPG -> HTL (2, intra-linker GPGPG, HTL_A
    # before HTL_B by nProm) -> inter-block GPGPG -> CTL (2, intra-linker
    # AAY, CTL_A before CTL_B because netcleaveMatch=True is prioritized
    # over raw promiscuity).
    EXPECTED_CONSTRUCT = (
        BCELL_GOOD + 'KK' + BCELL_ALLERGEN + 'GPGPG'
        + HTL_A + 'GPGPG' + HTL_B + 'GPGPG' + CTL_A + 'AAY' + CTL_B
    )

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls.protBcellROIs = cls._buildAnnotatedSet(
            seqName='BCELL_SEQ', windows=[(1, len(cls.BCELL_GOOD)), (len(cls.BCELL_GOOD) + 16, 60)],
            attrBuilder=lambda i: {
                '_algpredVerdict': String('Non-Allergen' if i == 0 else 'Allergen'),
                '_hasGlycoSequon': Boolean(False),
                '_glycoSequonSummary': String(''),
                '_meanScore': Float(0.8 if i == 0 else 0.2),
                '_maxScore': Float(0.9 if i == 0 else 0.3),
            },
            fullSeq=lambda: cls.BCELL_GOOD + 'X' * 15 + cls.BCELL_ALLERGEN,
        )

    @classmethod
    def _buildAnnotatedSet(cls, seqName, windows, attrBuilder, fullSeq):
        """Seeds a SetOfSequenceROIs via ProtImportSequence+ProtDefineSeqROI, then rebuilds
        it (full re-append, not update()) with the given dynamic attributes attached BEFORE
        append -- a Scipion Set fixes its column schema from the first appended item's
        attributes, so every item needs the exact same attribute set from the start."""
        seq = fullSeq()
        protImport = cls.newProtocol(
            ProtImportSequence, inputSequenceName=seqName, inputSequenceDescription=seqName,
            inputRawSequence=seq,
        )
        cls.proj.launchProtocol(protImport, wait=True)

        inROIs = '\n'.join(
            '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
                i, start, end, seq[start - 1:end]
            )
            for i, (start, end) in enumerate(windows, 1)
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(protImport)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')
        cls.proj.launchProtocol(protDefSeqROIs, wait=True)

        sqlitePath = protDefSeqROIs.outputROIs.getFileName()
        oldItems = [roi.clone() for roi in protDefSeqROIs.outputROIs]
        os.remove(sqlitePath)
        rebuilt = SetOfSequenceROIs(filename=sqlitePath)
        for i, roi in enumerate(oldItems):
            for attrName, value in attrBuilder(i).items():
                setattr(roi, attrName, value)
            rebuilt.append(roi)
        rebuilt.write()
        return protDefSeqROIs

    def test(self):
        # HTL and CTL fixtures built here (not in setUpClass) since their
        # window sequences embed the core9aa values directly as the FULL
        # window (a realistic simplification: in the real pipeline the
        # window is often wider than the core, but only core9aa is
        # inserted into the construct, so the window content beyond it is
        # irrelevant to this test).
        protHtlROIs = self._buildAnnotatedSet(
            seqName='HTL_SEQ', windows=[(1, 9), (20, 28)],
            attrBuilder=lambda i: {
                '_core9aa': String(self.HTL_A if i == 0 else self.HTL_B),
                '_nPromiscuousAlleles': Integer(5 if i == 0 else 3),
                '_nAllelesEvaluated': Integer(27),
                '_minRankEl': Float(0.5 if i == 0 else 1.0),
            },
            fullSeq=lambda: self.HTL_A + 'X' * 10 + self.HTL_B,
        )
        protCtlROIs = self._buildAnnotatedSet(
            seqName='CTL_SEQ', windows=[(1, 9), (20, 28)],
            attrBuilder=lambda i: {
                '_core9aa': String(self.CTL_A if i == 0 else self.CTL_B),
                '_nPromiscuousAlleles': Integer(4 if i == 0 else 6),
                '_nAllelesEvaluated': Integer(23),
                '_minRankEl': Float(0.8 if i == 0 else 0.3),
                '_netcleaveCTermMatch': Boolean(i == 0),
                '_netcleaveCTermScore': Float(0.69 if i == 0 else 0.0),
            },
            fullSeq=lambda: self.CTL_A + 'X' * 10 + self.CTL_B,
        )

        protAssembly = self.newProtocol(ProtEpitopeConstructAssembly)
        protAssembly.bcellROIs.set(self.protBcellROIs)
        protAssembly.bcellROIs.setExtended('outputROIs')
        protAssembly.htlROIs.set(protHtlROIs)
        protAssembly.htlROIs.setExtended('outputROIs')
        protAssembly.ctlROIs.set(protCtlROIs)
        protAssembly.ctlROIs.setExtended('outputROIs')
        self.launchProtocol(protAssembly, wait=True)

        outROIs = getattr(protAssembly, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), 1)

        roi = list(outROIs)[0]
        self.assertEqual(roi.getROISequence(), self.EXPECTED_CONSTRUCT)
        self.assertEqual(roi.getROIIdx(), 1)
        self.assertEqual(roi.getROIIdx2(), len(self.EXPECTED_CONSTRUCT))

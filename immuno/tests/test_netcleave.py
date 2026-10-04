import os

from pwchem.objects import SetOfSequenceROIs
from pwchem.protocols import ProtDefineSeqROI
from pwem.protocols import ProtImportSequence
from pyworkflow.object import String
from pyworkflow.tests import BaseTest, setupTestProject

from ..protocols import ProtNetCleavePrediction


class TestNetCleavePrediction(BaseTest):
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

    # Window 1 (562-584): real 23aa GP120 fragment. '_core9aa' is manually
    # attached below to its first 9 residues ('RAIEAQQHL', the real MHC-I
    # core found by the scipion-chem-netmhcpan test fixture), simulating
    # what an upstream MHC-I promiscuity protocol's output ROI looks like.
    # Expected c-term cleavage position = 0 (offset) + 9 (core len) + 1 = 10,
    # a REAL value confirmed by running the actual local NetCleave binary on
    # this exact window (score 0.694426, not estimated).
    #
    # Window 2 (745-765): no '_core9aa' attached -- exercises the fallback
    # (whole ROI treated as its own "core"), which structurally can never
    # match (the required cleavage position would fall past the window's
    # last tested position, confirmed for real): a genuine negative case.
    WINDOW_WITH_CORE = (562, 584)
    CORE_SEQ = 'RAIEAQQHL'
    EXPECTED_SCORE = 0.694426
    WINDOW_NO_CORE = (745, 765)

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

        cls.protSeedROIs = cls._runDefSeqROIs(cls.protImportSeq)
        cls._waitOutput(cls.protSeedROIs, 'outputROIs', sleepTime=5)

        cls._attachCore9aa()

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
        windows = [cls.WINDOW_WITH_CORE, cls.WINDOW_NO_CORE]
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

    @classmethod
    def _attachCore9aa(cls):
        """Simulates an upstream MHC-I promiscuity protocol's '_core9aa' attribute
        on the ROI covering WINDOW_WITH_CORE, so this test can exercise the
        real match path without needing to chain a full netmhcpan run.

        NOTE: Set.update() cannot add a new sqlite column to an
        already-persisted set (only Set.append() creates columns for new
        dynamic attributes, at append time) -- so the whole set is rebuilt
        from scratch (same file path, fresh rows) instead of patched in
        place.
        """
        sqlitePath = cls.protSeedROIs.outputROIs.getFileName()
        oldItems = [roi.clone() for roi in cls.protSeedROIs.outputROIs]

        # Delete the old file first: instantiating a Set with a filename
        # that already exists connects to it as-is (does not truncate it),
        # which would conflict with re-appending items that carry the same
        # pre-existing objIds.
        os.remove(sqlitePath)
        rebuilt = SetOfSequenceROIs(filename=sqlitePath)
        for roi in oldItems:
            # A Set fixes its column schema from the first appended item's
            # dynamic attributes: '_core9aa' must be set (even if blank) on
            # EVERY item, or later inserts crash with AttributeError.
            roi._core9aa = String(cls.CORE_SEQ if roi.getROIIdx() == cls.WINDOW_WITH_CORE[0] else '')
            rebuilt.append(roi)
        rebuilt.write()

    def test(self):
        protNetCleave = self.newProtocol(ProtNetCleavePrediction)
        protNetCleave.inputROIs.set(self.protSeedROIs)
        protNetCleave.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protNetCleave, wait=True)

        outROIs = getattr(protNetCleave, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), 2)

        # NOTE: iterating a Scipion SetOfXXX reuses the same underlying
        # Python object per row -- clone before storing for later use.
        roisByStart = {roi.getROIIdx(): roi.clone() for roi in outROIs}

        withCore = roisByStart[self.WINDOW_WITH_CORE[0]]
        self.assertTrue(withCore._netcleaveCTermMatch.get())
        self.assertAlmostEqual(withCore._netcleaveCTermScore.get(), self.EXPECTED_SCORE, places=5)

        noCore = roisByStart[self.WINDOW_NO_CORE[0]]
        self.assertFalse(noCore._netcleaveCTermMatch.get())

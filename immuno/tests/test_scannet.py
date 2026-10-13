from pwem.protocols import ProtImportPdb
from pyworkflow.tests import BaseTest, setupTestProject

from ..protocols import ProtScanNetPrediction

# 7c4s (antibody Fab + antigen complex, 6 chains) downloaded live from the
# RCSB PDB database at test time -- not a bundled/hardcoded local file, so
# the test is portable across machines.
_TEST_PDB_ID = '7c4s'


class TestScanNetPrediction(BaseTest):
    # Real reference values from a real local Docker run of ScanNet
    # (jertubiana/scannet image) against 7c4s.pdb (6 chains), followed by
    # this plugin's own adaptive-percentile (90th) + sliding-window (9aa,
    # max 2 below-threshold residues, min length 9) epitope mapping -- not
    # estimated. 7c4s has 2 pairs of near-identical chains (L/K, H/J), so
    # their regions are expected to match closely (real antibody Fab
    # symmetry), not a coincidence in this fixture.
    EXPECTED = sorted([
        ('L', 163, 172, 'WTDQDSKDST', 0.2722),
        ('H', 37, 47, 'VKQTPVHGLEW', 0.226455),
        ('A', 2, 26, 'QETLREHYQYVGKLAGRLKEASEGS', 0.15128),
        ('K', 38, 49, 'QKPEQSPKLLIY', 0.210333),
        ('K', 163, 172, 'WTDQDSKDST', 0.2712),
        ('J', 37, 47, 'VKQTPVHGLEW', 0.233455),
        ('B', 6, 25, 'REHYQYVGKLAGRLKEASEG', 0.1742),
    ])

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls.protImportPdb = cls._runImportPdb()

    @classmethod
    def _runImportPdb(cls):
        protImportPdb = cls.newProtocol(ProtImportPdb, inputPdbData=0, pdbId=_TEST_PDB_ID)
        cls.proj.launchProtocol(protImportPdb, wait=True)
        return protImportPdb

    def runScannet(self):
        protScanNet = self.newProtocol(ProtScanNetPrediction)
        protScanNet.inputStructure.set(self.protImportPdb)
        protScanNet.inputStructure.setExtended('outputPdb')
        self.launchProtocol(protScanNet, wait=True)
        return protScanNet

    def test(self):
        protScanNet = self.runScannet()

        outROIs = getattr(protScanNet, 'outputROIs', None)
        self.assertIsNotNone(outROIs)

        # NOTE: iterating a Scipion SetOfXXX reuses the same underlying
        # Python object per row -- read every needed value out immediately
        # per iteration instead of storing object references for later.
        got = sorted(
            (roi._sequence.getId().replace('chain_', ''), roi.getROIIdx(), roi.getROIIdx2(), roi.getROISequence(),
             round(roi._meanScore.get(), 6))
            for roi in outROIs
        )
        gotRounded = [(g[0], g[1], g[2], g[3]) for g in got]
        expectedRounded = [(e[0], e[1], e[2], e[3]) for e in self.EXPECTED]
        self.assertEqual(gotRounded, expectedRounded)
        for got_row, expected_row in zip(got, self.EXPECTED):
            self.assertAlmostEqual(got_row[4], expected_row[4], places=4)

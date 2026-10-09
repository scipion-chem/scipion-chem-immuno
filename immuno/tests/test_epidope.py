from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import ProtEpiDopePrediction


class TestEpiDopePrediction(BaseTest):
    """Requires a local EpiDope installation (see Plugin.addEpiDopePackage)."""

    NAME = 'CSP_NANP'
    DESCRIPTION = 'Plasmodium falciparum circumsporozoite protein, central NANP repeat region'
    AMINOACIDSSEQ = 'DPNANPNVDPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNA'

    # The NANP repeat is the immunodominant B-cell epitope of CSP, so every
    # residue of this fragment scores above EpiDope's default threshold and
    # the whole window collapses into a single region. Reference values
    # obtained by running the locally installed EpiDope CLI on this exact
    # fragment: per-residue scores 0.854-0.961.
    EXPECTED_ROIS = [(1, len(AMINOACIDSSEQ), AMINOACIDSSEQ)]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.DESCRIPTION,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    def runEpiDope(self):
        protEpiDope = self.newProtocol(ProtEpiDopePrediction)
        protEpiDope.inputSequence.set(self.protImportSeq)
        protEpiDope.inputSequence.setExtended('outputSequence')
        self.launchProtocol(protEpiDope, wait=True)
        return protEpiDope

    def test(self):
        protEpiDope = self.runEpiDope()

        outROIs = getattr(protEpiDope, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = [(roi.getROIIdx(), roi.getROIIdx2(), roi.getROISequence()) for roi in outROIs]
        self.assertEqual(got, self.EXPECTED_ROIS)

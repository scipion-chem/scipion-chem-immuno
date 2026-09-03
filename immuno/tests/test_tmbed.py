from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import ProtTMbedPredict


class TestTMbedPredict(BaseTest):
    NAME = 'VDAC1_P21796'
    DESCRIPTION = 'Voltage-dependent anion-selective channel protein 1, UniProt P21796'
    AMINOACIDSSEQ = (
        'MAVPPTYADLGKSARDVFTKGYGFGLIKLDLKTKSENGLEFTSSGSANTETTKVTGSLETKYRWTEYGLTFTEKWNTD'
        'NTLGTEITVEDQLARGLKLTFDSSFSPNTGKKNAKIKTGYKREHINLGCDMDFDIAGPSIRGALVLGYEGWLAGYQMN'
        'FETAKSRVTQSNFAVGYKTDEFQLHTNVNDGTEFGGSIYQKVNKKLETAVNLAWTAGNSNTRFGIAAKYQIDPDACFS'
        'AKVNNSSLIGLGYTQTLKPGIKLTLSALLDGKNVNAGGHKLGLGLEFQA'
    )

    # Real reference values: obtained by running the locally installed
    # TMbed CLI directly on this exact sequence (--out-format 1, CPU,
    # --model-dir pointing at the local ProtT5 weights), then collapsing the
    # per-residue class string with extract_masking_regions(). VDAC1 is a
    # 19-stranded beta-barrel outer mitochondrial membrane channel with no
    # signal peptide, so an all-'TM_beta_strand' result is the expected
    # biological shape, not an artifact.
    EXPECTED = sorted([
        (27, 32, 'TM_beta_strand'), (41, 46, 'TM_beta_strand'), (56, 62, 'TM_beta_strand'),
        (70, 74, 'TM_beta_strand'), (81, 88, 'TM_beta_strand'), (96, 102, 'TM_beta_strand'),
        (112, 119, 'TM_beta_strand'), (124, 131, 'TM_beta_strand'), (137, 145, 'TM_beta_strand'),
        (150, 156, 'TM_beta_strand'), (167, 173, 'TM_beta_strand'), (191, 196, 'TM_beta_strand'),
        (204, 209, 'TM_beta_strand'), (219, 225, 'TM_beta_strand'), (233, 237, 'TM_beta_strand'),
        (243, 248, 'TM_beta_strand'), (257, 262, 'TM_beta_strand'), (275, 280, 'TM_beta_strand'),
    ])

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

    def runTMbed(self):
        protTMbed = self.newProtocol(ProtTMbedPredict, useGpu=False, numThreads=4)
        protTMbed.inputSequence.set(self.protImportSeq)
        protTMbed.inputSequence.setExtended('outputSequence')
        self.launchProtocol(protTMbed, wait=True)
        return protTMbed

    def test(self):
        protTMbed = self.runTMbed()

        outROIs = getattr(protTMbed, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = sorted((roi.getROIIdx(), roi.getROIIdx2(), roi.getType()) for roi in outROIs)
        self.assertEqual(got, self.EXPECTED)

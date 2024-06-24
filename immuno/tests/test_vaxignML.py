# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import setupTestProject, DataSet, BaseTest

from pwchem.protocols import ProtChemImportSetOfSequences
from pwchem.utils import assertHandle

from ..protocols import ProtVaxignMLEpitopeEvaluation

class TestVaxignML(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.ds = DataSet.getDataSet('model_building_tutorial')

	@classmethod
	def _runImportSeqs(cls):
		protImportSeqs = cls.newProtocol(
			ProtChemImportSetOfSequences,
			multiple=True,
			filesPath=cls.ds.getFile('Sequences/'), filesPattern='*_A_mutated.fasta')
		cls.launchProtocol(protImportSeqs)
		cls.protImportSeqs = protImportSeqs
		return protImportSeqs

	def _runVaxignML(self, protSeqs):
		protSel = self.newProtocol(ProtVaxignMLEpitopeEvaluation,
															 inputSource=0)

		protSel.inputSequences.set(protSeqs)
		protSel.inputSequences.setExtended('outputSequences')

		self.proj.launchProtocol(protSel, wait=False)
		return protSel

	def test(self):
		protSeqs = self._runImportSeqs()
		self._waitOutput(protSeqs, 'outputSequences', sleepTime=10)
		protVax = self._runVaxignML(protSeqs)
		self._waitOutput(protVax, 'outputSequences', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protVax, 'outputSequences', None))

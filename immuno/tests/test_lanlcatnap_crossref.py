# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import unittest

from pwchem.tests import TestDefineMultiEpitope
from pwchem.utils import assertHandle

from immuno import Plugin as immunoPlugin
from ..constants import LANL_AB_ALL_PATH
from ..protocols import ProtLANLCATNAPCrossref

# Only this plugin's test is exported. TestDefineMultiEpitope is imported as a base
# class from pwchem and must not be collected as part of this plugin's suite.
__all__ = ['TestLANLCATNAPCrossref']


@unittest.skipUnless(
	os.path.isfile(immunoPlugin.getVar(LANL_AB_ALL_PATH) or ''),
	'LANL_AB_ALL_PATH is not configured. The LANL bnAb reference database is downloaded '
	'manually, see ProtLANLCATNAPCrossref and this plugin\'s README.')
class TestLANLCATNAPCrossref(TestDefineMultiEpitope):
	"""Requires LANL_AB_ALL_PATH to be configured, since the reference database is
	downloaded manually and cannot be redistributed. Skipped otherwise, see the class
	decorator."""

	@classmethod
	def _runLANLCATNAP(cls, protROIs):
		protCrossref = cls.newProtocol(ProtLANLCATNAPCrossref)
		protCrossref.inputROIs.set(protROIs)
		protCrossref.inputROIs.setExtended('outputROIs')

		cls.proj.launchProtocol(protCrossref, wait=False)
		return protCrossref

	def test(self):
		protsROIs = self._runDefSeqROIs(inProt=self.protImportSeq)
		self._waitOutput(protsROIs, 'outputROIs', sleepTime=5)

		protCrossref = self._runLANLCATNAP(protsROIs)
		self._waitOutput(protCrossref, 'outputROIs', sleepTime=5)
		assertHandle(self.assertIsNotNone, getattr(protCrossref, 'outputROIs', None),
								 cwd=protCrossref.getWorkingDir())

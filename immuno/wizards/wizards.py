# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
"""

from pwchem.wizards import AddElementWizard

from ..protocols import ProtIIITDEpitopeSelection, ProtIIITDEvaluations

class AddIIITDElement(AddElementWizard):
  """Add filter expression in ligand filter protocol"""
  _targets, _inputs, _outputs = [], {}, {}

  def buildSDicStr(self, sDic):
    return '\n'.join([f'{i+1}) "{sName}": {sd}' for i, (sName, sd) in enumerate(sDic.items())])

  def show(self, form, *params):
    _, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    prevSDic = protocol.parseElementsDic()
    selectorDic = protocol.buildElementDic()

    sKey, sValue = list(selectorDic.items())[0]
    if sKey not in prevSDic and sValue not in prevSDic.values():
      prevSDic.update(selectorDic)
      towrite = self.buildSDicStr(prevSDic)
      form.setVar(outputParam[0], towrite)
    else:
      print('The software name and values cannot be repeated')


AddIIITDElement().addTarget(protocol=ProtIIITDEpitopeSelection,
                             targets=['addSel'],
                             inputs=['inSels'],
                             outputs=['inSels'])

AddIIITDElement().addTarget(protocol=ProtIIITDEvaluations,
                             targets=['addEval'],
                             inputs=['inEvals'],
                             outputs=['inEvals'])

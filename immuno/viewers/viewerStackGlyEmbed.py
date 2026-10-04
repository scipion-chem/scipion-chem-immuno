# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
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
Viewer for ProtStackGlyEmbedPrediction: a per-site bar plot (region, verdict,
score) plus a filtered table view -- the raw SetOfSequenceROIs table buries
the glyco info among sequence/id columns.
"""

import re

from pyworkflow.protocol import params, Protocol
import pyworkflow.viewer as pwviewer
from pwem.viewers import EmPlotter, ObjectView, showj

from ..protocols import ProtStackGlyEmbedPrediction

# Matches one ';'-joined entry from ProtStackGlyEmbedPrediction.createOutputStep:
# '{pos}:{verdict}({score:.3f})' or '{pos}:{verdict}' (no score) when not evaluated.
_SITE_RE = re.compile(r'(\d+):([^(;]+?)(?:\(([\d.]+)\))?(?:;|$)')

VERDICT_COLORS = {
    'Glycosylated': 'tab:green',
    'Not glycosylated': 'tab:red',
    'Not evaluated': 'tab:gray',
}


class ViewerStackGlyEmbed(pwviewer.ProtocolViewer):
    """Visualizes StackGlyEmbed's per-sequon glycosylation calls: a bar plot
    (site position vs score, colored by verdict) and a filtered ROI table."""

    _label = 'StackGlyEmbed viewer'
    _targets = [ProtStackGlyEmbedPrediction]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def _defineParams(self, form):
        form.addSection(label='Visualization of StackGlyEmbed results')
        group = form.addGroup('Per-site glycosylation plot')
        group.addParam('displayGlycoPlot', params.LabelParam,
                        label='Display score by site: ',
                        help='Bar plot of every scanned N-X-[S/T] sequon: position on the parent '
                             'sequence (x axis), StackGlyEmbed score (y axis), colored by verdict '
                             '(green=Glycosylated, red=Not glycosylated, gray=Not evaluated).')

        group = form.addGroup('ROI table')
        group.addParam('displayROIsTable', params.LabelParam,
                        label='Display ROI table (region/state/score): ',
                        help='Output ROIs table, filtered to the region indices and glycosylation '
                             'columns (hides the sequence/id columns showj otherwise leads with).')

    def _getVisualizeDict(self):
        return {
            'displayGlycoPlot': self._showGlycoPlot,
            'displayROIsTable': self._showROIsTable,
        }

    # ---------------------------------- Helpers -----------------------------------

    def _getSites(self):
        """Flattens every ROI's ``_glycoSequonSummary`` into one row per sequon site.

        Returns:
            List of (roi_label, position, verdict, score_or_None) tuples,
            sorted by position. ``roi_label`` is ``'{roiIdx}-{roiIdx2}'``,
            kept for the table view (a site's absolute position already
            identifies it uniquely on the plot).
        """
        outROIs = getattr(self.getProtocol(), 'outputROIs', None)
        if outROIs is None:
            return []

        sites = []
        for roi in outROIs:
            roiLabel = f'{roi.getROIIdx()}-{roi.getROIIdx2()}'
            summary = roi._glycoSequonSummary.get() or ''
            for match in _SITE_RE.finditer(summary):
                pos, verdict, score = match.group(1), match.group(2), match.group(3)
                sites.append((roiLabel, int(pos), verdict, float(score) if score is not None else None))
        return sorted(sites, key=lambda row: row[1])

    def getProtocol(self):
        if hasattr(self, 'protocol') and isinstance(self.protocol, Protocol):
            return self.protocol

    # ---------------------------------- Views -----------------------------------

    def _showGlycoPlot(self, paramName=None):
        sites = self._getSites()
        if not sites:
            return []

        plotter = EmPlotter(x=1, y=1, windowTitle='StackGlyEmbed: score by site')
        ax = plotter.createSubPlot('N-glycosylation sequon scores', 'Sequence position (Asn)', 'Score')

        positions = [pos for _, pos, _, _ in sites]
        scores = [score if score is not None else 0.0 for _, _, _, score in sites]
        colors = [VERDICT_COLORS.get(verdict, 'tab:gray') for _, _, verdict, _ in sites]

        ax.bar([str(p) for p in positions], scores, color=colors)
        ax.axhline(0.5, color='black', linestyle='--', linewidth=0.8, label='0.5 threshold')
        ax.set_ylim(0, 1)
        ax.tick_params(axis='x', rotation=90)

        from matplotlib.patches import Patch
        handles = [Patch(color=c, label=v) for v, c in VERDICT_COLORS.items()]
        ax.legend(handles=handles, loc='upper right', fontsize=6)

        plotter.show()
        return [plotter]

    def _showROIsTable(self, paramName=None):
        outROIs = getattr(self.getProtocol(), 'outputROIs', None)
        if outROIs is None:
            return []

        viewParams = {
            showj.MODE: showj.MODE_TABLE,
            showj.ORDER: '_roiIdx,_roiIdx2,_hasGlycoSequon,_glycoSequonSummary',
            showj.VISIBLE: '_roiIdx,_roiIdx2,_hasGlycoSequon,_glycoSequonSummary',
        }
        return [ObjectView(self._project, outROIs.strId(), outROIs.getFileName(), viewParams=viewParams)]

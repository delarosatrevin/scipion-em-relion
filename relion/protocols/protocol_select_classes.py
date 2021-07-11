# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from emtable import Table
import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW
from pwem.protocols import ProtProcessParticles
from pwem.objects import SetOfClasses2D

from relion import Plugin
from ..constants import RELION_PYTHON


class ProtRelionSelectClasses2D(ProtProcessParticles):
    """
    Relion protocol to auto-select 2D class averages.
    """
    _label = '2D class ranker'
    _devStatus = NEW

    @classmethod
    def isDisabled(cls):
        return not Plugin.IS_GT31()

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {'cls_selection': self._getExtraPath('backup_selection.star')}
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputProtocol', params.PointerParam,
                      pointerClass='ProtRelionClassify2D',
                      label="Input Relion 2D classification",
                      important=True)
        form.addParam('minThreshold', params.FloatParam, default=0.5,
                      label='Min. threshold for auto-selection',
                      help='Only classes with a predicted threshold '
                           'above this value will be selected.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('runSelectStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def runSelectStep(self):
        inputProt = self.inputProtocol.get()
        inputProt._initialize()
        fnOptimiser = inputProt._getFileName('optimiser',
                                             iter=inputProt._lastIter())
        params = " --opt %s --o %s --min_score %s" % (fnOptimiser,
                                                      self._getExtraPath(),
                                                      self.minThreshold.get())
        params += " --fn_sel_parts particles.star"
        params += " --fn_sel_classavgs class_averages.star"
        params += " --fn_root rank --do_granularity_features"
        params += " --auto_select"
        params += " --python %s" % Plugin.getVar(RELION_PYTHON)

        self.runJob("relion_class_ranker", params)

    def createOutputStep(self):
        table = Table(fileName=self._getFileName('cls_selection'))
        self._clsSelection = table.getColumnValues('rlnSelected')

        inputClasses = self.inputProtocol.get().outputClasses
        output = SetOfClasses2D.create(self._getExtraPath())
        output.copyInfo(inputClasses)
        output.appendFromClasses(inputClasses, filterClassFunc=self._appendClass)

        self._defineOutputs(outputClasses=output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _validate(self):
        errors = []
        if Plugin.getVar(RELION_PYTHON) is None:
            errors.append("%s is not defined in the Scipion configuration!\n"
                          "Please set it to point to a Python that "
                          "includes torch and numpy modules." % RELION_PYTHON)

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _appendClass(self, item):
        return False if not self._clsSelection[item.getObjId()-1] else True

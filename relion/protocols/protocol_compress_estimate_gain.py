# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# ******************************************************************************

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.protocol import STEPS_SERIAL
from pyworkflow.constants import PROD
from pwem.protocols import ProtProcessMovies

import relion
import relion.convert as convert


class ProtRelionCompressEstimateGain(ProtProcessMovies):
    """
    Using *relion_convert_to_tiff* to estimate the gain that can be used
    for better compression.
    """
    _label = 'estimate gain to compress'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtProcessMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=pwutils.Message.LABEL_INPUT)

        form.addParam('inputMovies', params.PointerParam,
                      pointerClass='SetOfMovies',
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of movies to be used in the gain '
                           'estimation.')

        form.addParam('moviesSubset', params.IntParam, default=0,
                      label='Subset',
                      help="Use a subset of the movies for the estimation. "
                           "If 0, all input movies will be used. ")

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputMovies.getObjId(),
                                 self.moviesSubset.get())

        self._insertFunctionStep('estimateGainStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self, moviesId, subset):
        self.info("Relion version:")
        self.runJob("relion_convert_to_tiff --version", "", numberOfMpi=1)
        self.info("Detected version from config: %s"
                  % relion.Plugin.getActiveVersion())

        moviesList = self.inputMovies.get()
        if subset > 0:
            moviesList = [m.clone()
                          for i, m in enumerate(moviesList) if i < subset]

        micStar = self._getTmpPath('input_movies.star')
        tmp = self._getTmpPath()
        writer = convert.createWriter(rootDir=tmp, outputDir=tmp)
        writer.writeSetOfMovies(moviesList, micStar)

    def estimateGainStep(self):
        args = "--i input_movies.star --o ../ "
        args += "--estimate_gain --j %s " % self.numberOfThreads

        self.runJob('relion_convert_to_tiff', args, cwd=self._getTmpPath())

    def createOutputStep(self):
        pass

    def _stepsCheck(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        errors = []

        return errors

    def getOutputGain(self):
        return self._getPath('gain_estimate.bin')

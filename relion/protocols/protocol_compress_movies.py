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

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pwem.protocols import ProtAlignMovies
from pwem.objects import MovieAlignment
from pyworkflow.protocol import STEPS_PARALLEL

import relion


class ProtRelionCompressMovies(ProtAlignMovies):
    """
    Using *relion_convert_to_tiff* to compress a set of movies.
    """
    _label = 'compress movies'
    _devStatus = PROD

    OP_COMPRESS = 0
    OP_ESTIMATE = 1

    def __init__(self, **kwargs):

        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):

        form.addParam('inputGainProt', params.PointerParam, allowsNull=True,
                      pointerClass='ProtRelionCompressEstimateGain',
                      label='Input gain (Optional)',
                      help='Provide as input a compress estimate protocol '
                           'from where the estimated gain file will be taken.')

        group = form.addGroup("TIFF Options")

        group.addParam('compression', params.EnumParam, default=1,
                       choices=['none', 'auto', 'zip', 'lzw'],
                       label='Compression type')

        group.addParam('deflateLevel', params.IntParam, default=6,
                       label='Deflate level',
                       help="deflate level. 1 (fast) "
                            "to 9 (slowest but best compression)")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions -------------------------------
    def _convertInputStep(self):
        self.info("Relion version:")
        self.runJob("relion_convert_to_tiff --version", "", numberOfMpi=1)
        self.info("Detected version from config: %s"
                  % relion.Plugin.getActiveVersion())
        # Create a local link to the input gain file if necessary
        inputGain = self.getInputGain()
        if inputGain:
            tmpGain = self._getTmpPath('gain_estimate.bin')
            pwutils.createLink(inputGain, tmpGain)
            pwutils.createLink(inputGain.replace('.bin', '_reliablity.bin'),
                               tmpGain.replace('.bin', '_reliablity.bin'))
        ProtAlignMovies._convertInputStep(self)

    def _processMovie(self, movie):
        fn = movie.getAttributeValue('_originalFileName', movie.getFileName())
        baseName = os.path.basename(fn)
        compression = self.getEnumText('compression')
        pwutils.createLink(fn, self._getTmpPath(baseName))
        args = "--i %s --o ../extra/ " % baseName
        args += "--compression %s " % compression
        # TODO: Check if deflateLevel is only valid for zip (deflate)
        if compression == 'zip':  # deflate
            args += "--deflate_level %d" % self.deflateLevel

        if self.getInputGain():
            args += "--gain gain_estimate.bin "

        self.runJob('relion_convert_to_tiff', args, cwd=self._getTmpPath())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        # Check base validation before the specific ones for Motioncor
        errors = ProtAlignMovies._validate(self)

        return errors

    # ------------------------ Extra BASE functions ---------------------------
    def _getRelPath(self, baseName, refPath):
        return os.path.relpath(self._getExtraPath(baseName), refPath)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _createOutputMovies(self):
        return True

    def _createOutputMicrographs(self):
        return False

    def _createOutputWeightedMicrographs(self):
        return False

    def _savePsSum(self):
        return self.getAttributeValue('savePSsum', False)

    def _createOutputMovie(self, movie):
        n = movie.getNumberOfFrames()
        newMovie = movie.clone()
        movieName = pwutils.replaceBaseExt(movie.getFileName(), 'tif')
        newMovie.setFileName(self._getExtraPath(movieName))
        newMovie.setAlignment(MovieAlignment(first=1, last=n,
                                             xshifts=[0]*n, yshifts=[0]*n))
        return newMovie

    def getInputGain(self):
        inputGainProt = self.inputGainProt.get()
        if inputGainProt is not None:
            inputGain = inputGainProt.getOutputGain()
            if os.path.exists(inputGain):
                return inputGain
        return None

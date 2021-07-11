# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     Grigory Sharov     (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
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
from math import ceil
import json
import emtable as md

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.constants import PROD
import pyworkflow.utils as pwutils
from pwem.protocols import ProtAlignMovies
from pwem.objects import Image
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_SERIAL

import relion
from relion import Plugin
import relion.convert as convert
from relion.convert.convert31 import OpticsGroups
from .protocol_base import ProtRelionBase


class ProtRelionMotioncor(ProtAlignMovies, ProtRelionBase):
    """ Wrapper for the Relion's implementation of motioncor algorithm. """

    _label = 'motion correction'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        self.updatedSets = []
        self.isEER = False

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif', '.eer', '.gain'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        line = form.addLine('Frames for corrected SUM',
                            help='First and last frames to use in corrected '
                                 'average (starts counting at 1 and 0 as last '
                                 'means until the last frame in the movie). '
                                 'When using EER, the frames are not hardware '
                                 'frames, but fractions.')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')

        form.addParam('doDW', params.BooleanParam, default=True,
                      label='Do dose-weighting?',
                      help='If set to Yes, the averaged micrographs will be '
                           'dose-weighted.')

        form.addParam('saveNonDW', params.BooleanParam, default=False,
                      condition='doDW',
                      label='Save non-dose weighted as well?',
                      help='Aligned but non-dose weighted images are '
                           'sometimes useful in CTF estimation, although '
                           'there is no difference in most cases. Whichever '
                           'the choice, CTF refinement job is always done on '
                           'dose-weighted particles.')

        form.addParam('savePSsum', params.BooleanParam, default=False,
                      label='Save sum of power spectra?',
                      help='Sum of non-dose weighted power spectra '
                           'provides better signal for CTF estimation. '
                           'The power spectra can be used by CTFFIND4 '
                           'but not by GCTF.')
        form.addParam('dosePSsum', params.FloatParam, default=4.0,
                      condition='savePSsum',
                      label='Sum power spectra every e/A2',
                      help='McMullan et al. (Ultramicroscopy, 2015) '
                           'suggests summing power spectra every '
                           '4.0 e/A2 gives optimal Thon rings.')

        if Plugin.IS_GT31():
            form.addParam('saveFloat16', params.BooleanParam, default=False,
                          label="Write output in float16?",
                          expertLevel=params.LEVEL_ADVANCED,
                          lavel="Write output in float16?",
                          help="Relion can write output images in float16 "
                               "MRC (mode 12) format to save disk space. "
                               "By default, float32 format is used.")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      label="Compute PSD?",
                      help="If Yes, the protocol will compute for each "
                           "aligned micrograph the PSD using EMAN2.")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail with EMAN2 and keep it with the '
                           'micrograph object for visualization purposes.')

        form.addParam('extraParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="Extra parameters for Relion motion correction. ")

        form.addSection("Motion")
        form.addParam('bfactor', params.IntParam, default=150,
                      label='Bfactor',
                      help="The B-factor that will be applied to the "
                           "micrographs.")

        line = form.addLine('Number of patches',
                            help='Number of patches (in X and Y direction) to '
                                 'apply motion correction. If <= 2 then '
                                 'only global correction will be done.')
        line.addParam('patchX', params.IntParam, default=1, label='X')
        line.addParam('patchY', params.IntParam, default=1, label='Y')

        form.addParam('groupFrames', params.IntParam, default=1,
                      label='Group frames',
                      help="Average together this many frames before "
                           "calculating the beam-induced shifts.")

        form.addParam('binFactor', params.FloatParam, default=1.,
                      label='Binning factor',
                      help='Bin the micrographs this much by a windowing '
                           'operation in the Fourier Transform. Binning at '
                           'this level is hard to un-do later on, but may be '
                           'useful to down-scale super-resolution images. '
                           'Float-values may be used. Do make sure though '
                           'that the resulting micrograph size is even.')

        form.addParam('gainRot', params.EnumParam, default=0,
                      choices=['No rotation (0)',
                               ' 90 degrees (1)',
                               '180 degrees (2)',
                               '270 degrees (3)'],
                      label='Gain rotation',
                      help="Rotate the gain reference by this number times 90 "
                           "degrees clockwise in relion_display. This is the "
                           "same as -RotGain in MotionCor2. \n"
                           "Note that MotionCor2 uses a different convention "
                           "for rotation so it says 'counter-clockwise'.")

        form.addParam('gainFlip', params.EnumParam, default=0,
                      choices=['No flipping        (0)',
                               'Flip upside down   (1)',
                               'Flip left to right (2)'],
                      label='Gain flip',
                      help="Flip the gain reference after rotation. "
                           "This is the same as -FlipGain in MotionCor2. "
                           "0 means do nothing, 1 means flip Y (upside down) "
                           "and 2 means flip X (left to right).")

        form.addParam('defectFile', params.FileParam, allowsNull=True,
                      label='Defects file',
                      help='Location of a UCSF MotionCor2-style '
                           'defect text file or a defect map that '
                           'describe the defect pixels on the detector. '
                           'Each line of a defect text file should contain '
                           'four numbers specifying x, y, width and height '
                           'of a defect region. A defect map is an image '
                           '(MRC or TIFF), where 0 means good and 1 means '
                           'bad pixels. The coordinate system is the same '
                           'as the input movie before application of '
                           'binning, rotation and/or flipping.\n\n'
                           '_Note that the format of the defect text is '
                           'DIFFERENT from the defect text produced '
                           'by SerialEM!_\n One can convert a SerialEM-style '
                           'defect file into a defect map using IMOD '
                           'utilities e.g.:\n'
                           '*clip defect -D defect.txt -f tif movie.tif defect_map.tif*\n'
                           'See explanations in the SerialEM manual.\n'
                           'Leave empty if you do not have any defects, '
                           'or do not want to correct for defects on your detector.')

        form.addSection("EER")
        form.addParam('eerGroup', params.IntParam, default=32,
                      label='EER fractionation',
                      help="The number of hardware frames to group into one "
                           "fraction. This option is relevant only for Falcon4 "
                           "movies in the EER format. Falcon 4 operates at "
                           "248 frames/s.\nFractionate such that each fraction "
                           "has about 0.5 to 1.25 e/A2.")
        form.addParam('eerSampling', params.EnumParam, default=0,
                      expertLevel=cons.LEVEL_ADVANCED,
                      choices=['1', '2'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='EER upsampling',
                      help="EER upsampling (1 = 4K or 2 = 8K). 8K rendering is not "
                           "recommended by Relion. See "
                           "https://relion.readthedocs.io/en/latest/Reference/MovieCompression.html")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions -------------------------------
    def _convertInputStep(self):
        self.info("Relion version:")
        self.runJob("relion_run_motioncorr --version", "", numberOfMpi=1)
        self.info("Detected version from config: %s" % relion.Plugin.getActiveVersion())

        ProtAlignMovies._convertInputStep(self)

    def _processMovie(self, movie):
        movieFolder = self._getOutputMovieFolder(movie)
        inputStar = os.path.join(movieFolder,
                                 '%s_input.star' % self._getMovieRoot(movie))
        pwutils.makePath(os.path.join(movieFolder, 'output'))

        og = OpticsGroups.fromImages(self.inputMovies.get())
        writer = convert.createWriter(optics=og)
        # Let's use only the basename, since we will launch the command
        # from the movieFolder
        movie.setFileName(os.path.basename(movie.getFileName()))
        writer.writeSetOfMovies([movie], inputStar)

        # The program will run in the movie folder, so let's put
        # the input files relative to that
        args = "--i %s --o output/ " % os.path.basename(inputStar)
        args += "--use_own "
        args += "--first_frame_sum %d --last_frame_sum %d " % (self._getFrameRange())
        args += "--bin_factor %f --bfactor %d " % (self.binFactor, self.bfactor)
        args += "--angpix %0.5f " % (movie.getSamplingRate())
        args += "--patch_x %d --patch_y %d " % (self.patchX, self.patchY)
        args += "--group_frames %d " % self.groupFrames
        args += "--j %d " % self.numberOfThreads

        inputMovies = self.inputMovies.get()
        if inputMovies.getGain():
            args += '--gainref "%s" ' % inputMovies.getGain()
            args += '--gain_rot %d ' % self.gainRot
            args += '--gain_flip %d ' % self.gainFlip

        if self.defectFile.get():
            args += '--defect_file "%s" ' % self.defectFile.get()

        if self._savePsSum():
            args += '--grouping_for_ps %d ' % self._calcPsDose()

        if self.doDW:
            args += "--dose_weighting "
            preExp, dose = self._getCorrectedDose(inputMovies)
            # when using EER, the hardware frames are grouped
            if self.isEER:
                dose *= self.eerGroup.get()
            args += "--dose_per_frame %f " % dose
            args += "--preexposure %f " % preExp

            if self.saveNonDW:
                args += "--save_noDW "

        if self.isEER:
            args += "--eer_grouping %d " % self.eerGroup
            args += "--eer_upsampling %d " % (self.eerSampling.get() + 1)

        if Plugin.IS_GT31() and self.saveFloat16:
            args += "--float16 "

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        try:
            self._runProgram('relion_run_motioncorr', args, cwd=movieFolder)
            try:
                self._saveAlignmentPlots(movie, inputMovies.getSamplingRate())
                self._computeExtra(movie)
            except:
                self.error("ERROR: Extra work "
                           "(i.e plots, PSD, thumbnail) has failed for %s\n"
                           % movie.getFileName())

            self._moveFiles(movie)
        except:
            print("ERROR processing movie: ", movie.getFileName())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        errors = []
        inputMovies = self.inputMovies.get()

        # check if the first movie exists
        firstMovie = self.inputMovies.get().getFirstItem()
        if not os.path.exists(firstMovie.getFileName()):
            errors.append("The input movie files do not exist!!! "
                          "Since usually input movie files are symbolic links, "
                          "please check that links are not broken if you "
                          "moved the project folder. ")

        # check frames range
        _, lastFrame, _ = inputMovies.getFramesRange()

        #FIXME: Remove once EER support is enabled
        if pwutils.getExt(firstMovie.getFileName()) == ".eer":
            errors.append("EER support is not available yet in this plugin")

        #self.isEER = pwutils.getExt(firstMovie.getFileName()) == ".eer"
        if self.isEER:
            lastFrame //= self.eerGroup.get()

        msg = "Frames range to sum must be within %d - %d" % (1, lastFrame)
        if self.sumFrameN.get() == 0:
            self.sumFrameN.set(lastFrame)

        if not (1 <= self.sumFrame0 < lastFrame):
            errors.append(msg)
        elif not (self.sumFrameN <= lastFrame):
            errors.append(msg)
        elif not (self.sumFrame0 < self.sumFrameN):
            errors.append(msg)

        # check dose for DW
        acq = inputMovies.getAcquisition()
        if self.doDW:
            dose = acq.getDosePerFrame()
            if dose is None or dose < 0.00001:
                errors.append("Input movies do not contain the dose per frame, "
                              "dose-weighting can not be performed.")

        # check grouping in EER case
        if self.isEER and self.groupFrames.get() != 1:
            errors.append("*Group frames* option should be set to 1 when processing "
                          "EER movies. Please use *EER fractionation* option "
                          "instead.")

        # check eman2 plugin
        if self.doComputeMicThumbnail or self.doComputePSD:
            try:
                from pwem import Domain
                eman2 = Domain.importFromPlugin('eman2', doRaise=True)
            except:
                errors.append("EMAN2 plugin not found!\nComputing thumbnails "
                              "or PSD requires EMAN2 plugin and binaries installed.")

        if self.hasAttribute('saveFloat16') and self.saveFloat16:
            errors.append("MRC float16 format is not yet supported by XMIPP, "
                          "so you cannot use this option.")

        return errors

    # ------------------------ Extra BASE functions ---------------------------
    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _createOutputMovies(self):
        return True

    def _createOutputMicrographs(self):
        return not bool(self.doDW) or bool(self.saveNonDW)

    def _createOutputWeightedMicrographs(self):
        return bool(self.doDW)

    def _savePsSum(self):
        return self.getAttributeValue('savePSsum', False)

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)
        self._setMotionValues(movie, mic)
        if self._savePsSum():
            outPs = self._getExtraPath(self._getOutputMicPsName(movie))
            mic._powerSpectra = Image(location=outPs)
            mic._powerSpectra.setSamplingRate(self._calcPSSampling())

    def _setMotionValues(self, movie, mic):
        """ Parse motion values from the 'corrected_micrographs.star' file
        generated for each movie. """
        fn = self._getMovieExtraFn(movie, 'corrected_micrographs.star')
        table = md.Table(fileName=fn, tableName='micrographs')
        row = table[0]
        mic._rlnAccumMotionTotal = pwobj.Float(row.rlnAccumMotionTotal)
        mic._rlnAccumMotionEarly = pwobj.Float(row.rlnAccumMotionEarly)
        mic._rlnAccumMotionLate = pwobj.Float(row.rlnAccumMotionLate)

    def _getMovieShifts(self, movie, outStarFn=None):
        outStar = outStarFn or self._getMovieExtraFn(movie, '.star')
        first, last = self._getFrameRange()
        n = last - first + 1
        table = md.Table(fileName=outStar, tableName='global_shift')
        xShifts, yShifts = [], []

        for row in table:
            # Shifts are in pixels of the original (unbinned) movies
            xShifts.append(float(row.rlnMicrographShiftX))
            yShifts.append(float(row.rlnMicrographShiftY))
            if len(xShifts) == n:
                break

        return xShifts, yShifts

    def _getBinFactor(self):
        if not self.isEER:
            return self.binFactor.get()
        else:
            return self.binFactor.get() / (self.eerSampling.get() + 1)

    # --------------------------- UTILS functions -----------------------------
    def _getMovieOutFn(self, movie, suffix):
        movieBase = pwutils.removeBaseExt(movie.getFileName()).replace('.', '_')
        return os.path.join(self._getOutputMovieFolder(movie), 'output',
                            '%s%s' % (movieBase, suffix))

    def _getMovieExtraFn(self, movie, suffix):
        """ Return filenames in the extra directory with the prefix of this movie.
        Used to keep files associated with each micrograph.
        """
        movieBase = pwutils.removeBaseExt(movie.getFileName())
        return self._getExtraPath('%s%s' % (movieBase, suffix))

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd', 'png', extra=True)

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = Image(location=self._getPsdCorr(movie))
        if self.doComputeMicThumbnail:
            mic.thumbnail = Image(location=self._getOutputMicThumbnail(movie))

    def _computeExtra(self, movie):
        """ Compute thumbnail and PSD. """
        outMicFn = self._getMovieOutFn(movie, '.mrc')

        if self.doComputeMicThumbnail:
            self.computeThumbnail(outMicFn,
                                  outputFn=self._getOutputMicThumbnail(movie))
        if self.doComputePSD:
            self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))

    def _computePSD(self, inputFn, outputFn, scaleFactor=6):
        """ Generate a thumbnail of the PSD with EMAN2"""
        args = "%s %s " % (inputFn, outputFn)
        args += "--process=math.realtofft --meanshrink %s " % scaleFactor
        args += "--fixintscaling=sane"

        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2')
        from pyworkflow.utils.process import runJob
        runJob(self._log, eman2.Plugin.getProgram('e2proc2d.py'), args,
               env=eman2.Plugin.getEnviron())

        return outputFn

    def _moveFiles(self, movie):
        # It really annoying that Relion default names changes if you use DW or not
        # if use DW, the default name are DW and the others noDW
        if self.doDW:
            pwutils.moveFile(self._getMovieOutFn(movie, '.mrc'),
                             self._getExtraPath(self._getOutputMicWtName(movie)))
            if self.saveNonDW:
                pwutils.moveFile(self._getMovieOutFn(movie, '_noDW.mrc'),
                                 self._getExtraPath(self._getOutputMicName(movie)))
        else:
            pwutils.moveFile(self._getMovieOutFn(movie, '.mrc'),
                             self._getExtraPath(self._getOutputMicName(movie)))

        if self._savePsSum():
            pwutils.moveFile(self._getMovieOutFn(movie, '_PS.mrc'),
                             self._getExtraPath(self._getOutputMicPsName(movie)))

        # Keep some local files of this movie in the extra folder
        for suffix in ['.star', '.log']:
            pwutils.moveFile(self._getMovieOutFn(movie, suffix),
                             self._getMovieExtraFn(movie, suffix))

        suffix = 'corrected_micrographs.star'
        fn = os.path.join(self._getOutputMovieFolder(movie), 'output', suffix)
        pwutils.moveFile(fn, self._getMovieExtraFn(movie, suffix))

    def _saveAlignmentPlots(self, movie, pixSize):
        # Create plots and save as an image
        shiftsX, shiftsY = self._getMovieShifts(movie, self._getMovieOutFn(movie, '.star'))
        first, _ = self._getFrameRange()
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _createOutputMovie(self, movie):
        """ Overwrite this function to store the Relion's specific
        Motion model coefficients and Hot pixels.
        """
        m = ProtAlignMovies._createOutputMovie(self, movie)
        # Load local motion values only if the patches are more than one
        if self.patchX.get() > 2 and self.patchY.get() > 2:
            try:
                table = md.Table(fileName=self._getMovieExtraFn(movie, '.star'),
                                 tableName='local_motion_model')
                coeffs = [row.rlnMotionModelCoeff for row in table]
            except:
                print("Failed to parse local motion from: %s" %
                      os.path.abspath(self._getMovieExtraFn(movie, '.star')))
                coeffs = []  # Failed to parse the local motion
            m._rlnMotionModelCoeff = pwobj.String(json.dumps(coeffs))

        try:
            table = md.Table(fileName=self._getMovieExtraFn(movie, '.star'),
                             tableName='hot_pixels')
            hotPixels = [(row.rlnCoordinateX, row.rlnCoordinateY) for row in table]
        except:
            hotPixels = []
        m._rlnHotPixels = pwobj.String(json.dumps(hotPixels))

        return m

    def createOutputStep(self):
        # This method is re-implemented here because a bug in the base protocol
        # where the outputMicrographs is used without check if it is produced.
        # validate that we have some output movies
        if self._createOutputMovies():
            output = self.outputMovies
        elif self._createOutputMicrographs():
            output = self.outputMicrographs
        elif self._createOutputWeightedMicrographs():
            output = self.outputMicrographsDoseWeighted
        else:
            raise Exception("It does not seem like any output is produced!")

        inputSize = len(self.listOfMovies)
        outputSize = output.getSize()

        if outputSize == 0 and inputSize != 0:
            raise Exception("All movies failed, didn't create outputMicrographs."
                            "Please review movie processing steps above.")

        if outputSize < inputSize:
            self.warning(pwutils.yellowStr("WARNING - Failed to align %d movies."
                                           % (inputSize - outputSize)))

    def _updateOutputSet(self, outputName, outputSet,
                         state=pwobj.Set.STREAM_OPEN):
        """ Redefine this method to update optics info. """

        ogDict = {'rlnMicrographOriginalPixelSize': self.inputMovies.get().getSamplingRate()}

        if self.isEER:
            ogDict.update({'rlnEERUpsampling': self.eerSampling.get() + 1,
                           'rlnEERGrouping': self.eerGroup.get()})

        if outputName not in self.updatedSets:
            og = OpticsGroups.fromImages(outputSet)
            og.updateAll(**ogDict)
            og.toImages(outputSet)
            self.updatedSets.append(outputName)

        ProtAlignMovies._updateOutputSet(self, outputName, outputSet,
                                         state=state)

    def _calcPsDose(self):
        _, dose = self._getCorrectedDose(self.inputMovies.get())
        # when using EER, the hardware frames are grouped
        if self.isEER:
            dose *= self.eerGroup.get()
        dose_for_ps = round(self.dosePSsum.get() / dose)

        return 1 if dose_for_ps == 0 else dose_for_ps

    def _calcPSSampling(self):
        """ Adapted from relion 3.1 code. """
        movieSet = self.inputMovies.get()
        target_pixel_size = 1.4  # from CTFFIND 4.1
        ps = movieSet.getSamplingRate()
        ps_angpix = ps * self.binFactor.get()
        x, y, _ = movieSet.getDimensions()
        ps_size_square = min(x, y)
        if ps_angpix < target_pixel_size:
            nx_needed = ceil(ps_size_square * ps_angpix / target_pixel_size)
            nx_needed += nx_needed % 2
            ps_angpix = ps_size_square * ps_angpix / nx_needed
        return ps_angpix

    def _getOutputMicPsName(self, movie):
        """ Returns the name of the output PS
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic_PS.mrc'

    def _getFrameRange(self, n=None, prefix=None):
        # Reimplement this method
        return self.sumFrame0.get(), self.sumFrameN.get()


def createGlobalAlignmentPlot(meanX, meanY, first, pixSize):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []

    def px_to_ang(ax_px):
        y1, y2 = ax_px.get_ylim()
        x1, x2 = ax_px.get_xlim()
        ax_ang2.set_ylim(y1 * pixSize, y2 * pixSize)
        ax_ang.set_xlim(x1 * pixSize, x2 * pixSize)
        ax_ang.figure.canvas.draw()
        ax_ang2.figure.canvas.draw()

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax_px = figure.add_subplot(111)
    ax_px.grid()
    ax_px.set_xlabel('Shift x (px)')
    ax_px.set_ylabel('Shift y (px)')

    ax_ang = ax_px.twiny()
    ax_ang.set_xlabel('Shift x (A)')
    ax_ang2 = ax_px.twinx()
    ax_ang2.set_ylabel('Shift y (A)')

    i = first
    # The output _rlnMicrographShiftX/Y shifts relative to the first frame.
    # Unit is pixels of the original (unbinned) movies (Takanori, 2018)
    skipLabels = ceil(len(meanX) / 10.0)
    labelTick = 1

    for x, y in zip(meanX, meanY):
        sumMeanX.append(x)
        sumMeanY.append(y)
        if labelTick == 1:
            ax_px.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    # automatically update lim of ax_ang when lim of ax_px changes.
    ax_px.callbacks.connect("ylim_changed", px_to_ang)
    ax_px.callbacks.connect("xlim_changed", px_to_ang)

    ax_px.plot(sumMeanX, sumMeanY, color='b')
    ax_px.plot(sumMeanX, sumMeanY, 'yo')
    ax_px.plot(sumMeanX[0], sumMeanY[0], 'ro', markersize=10, linewidth=0.5)
    ax_px.set_title('Global frame alignment')

    plotter.tightLayout()

    return plotter

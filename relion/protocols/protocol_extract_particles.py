# **************************************************************************
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
# **************************************************************************

import os

from pyworkflow.object import Set, Integer
import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STATUS_FINISHED
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD

from pwem.protocols import ProtExtractParticles
from pwem.emlib.image import ImageHandler
from pwem.objects import Particle, Acquisition

import relion.convert
from relion import Plugin
from relion.convert.convert31 import OpticsGroups
from relion.constants import OTHER
from .protocol_base import ProtRelionBase


class ProtRelionExtractParticles(ProtExtractParticles, ProtRelionBase):
    """ Protocol to extract particles using a set of coordinates. """

    _label = 'particles extraction'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtExtractParticles.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _definePreprocessParams(self, form):
        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='This is size of the boxed particles (in pixels).')

        form.addParam('doRescale', params.BooleanParam, default=False,
                      label='Rescale particles?',
                      help='If set to Yes, particles will be re-scaled. '
                           'Note that the re-scaled size below will be in '
                           'the down-scaled images.')

        form.addParam('rescaledSize', params.IntParam, default=128,
                      validators=[params.Positive],
                      condition='doRescale',
                      label='Re-scaled size (px)',
                      help='Final size in pixels of the extracted particles. '
                           'The provided value should be an even number. ')

        if relion.Plugin.IS_GT31():
            form.addParam('saveFloat16', params.BooleanParam, default=False,
                          expertLevel=params.LEVEL_ADVANCED,
                          label="Write output in float16?",
                          help="Relion can write output images in float16 "
                               "MRC (mode 12) format to save disk space. "
                               "By default, float32 format is used.")

        form.addSection(label='Preprocess')

        form.addParam('doInvert', params.BooleanParam, default=None,
                      label='Invert contrast?',
                      help='Invert the contrast if your particles are black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize particles?',
                      help='If set to Yes, particles will be normalized in '
                           'the way RELION prefers it.')

        form.addParam('backDiameter', params.IntParam, default=-1,
                      condition='doNormalize',
                      label='Diameter background circle before scaling (px)',
                      help='Particles will be normalized to a mean value of '
                           'zero and a standard-deviation of one for all '
                           'pixels in the background area. The background area '
                           'is defined as all pixels outside a circle with '
                           'this given diameter in pixels (before rescaling). '
                           'When specifying a negative value, a default value '
                           'of 75% of the Particle box size will be used.')

        form.addParam('stddevWhiteDust', params.FloatParam, default=-1,
                      condition='doNormalize',
                      label='Stddev for white dust removal: ',
                      help='Remove very white pixels from the extracted '
                           'particles. Pixels values higher than this many '
                           'times the image stddev will be replaced with '
                           'values from a Gaussian distribution. \n'
                           'Use negative value to switch off dust removal.')

        form.addParam('stddevBlackDust', params.FloatParam, default=-1,
                      condition='doNormalize',
                      label='Stddev for black dust removal: ',
                      help='Remove very black pixels from the extracted '
                           'particles. Pixels values higher than this many '
                           'times the image stddev will be replaced with '
                           'values from a Gaussian distribution. \n'
                           'Use negative value to switch off dust removal.')

        self._defineStreamingParams(form)

        form.addParallelSection(threads=0, mpi=4)

    # -------------------------- INSERT steps functions -----------------------
    def _insertInitialSteps(self):
        self._setupBasicProperties()

        # used to convert micrographs if not in .mrc format
        self._ih = ImageHandler()

        # dict to map between micrographs and its coordinates file
        self._micCoordStarDict = {}

        # When no streaming, it doesn't make sense the default value of
        # batch size = 1, so let's use 0 to extract all micrographs at once
        if not self._isStreamOpen() and self._getStreamingBatchSize() == 1:
            self.info("WARNING: The batch size of 1 does not make sense when "
                      "not in streaming...changed value to 0 (extract all).")
            self.streamingBatchSize.set(0)

        return [self._insertFunctionStep('convertInputStep',
                                         self.getInputMicrographs().getObjId())]

    def _doNothing(self, *args):
        pass  # used to avoid some streaming functions

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self, micsId):
        self.info("Relion version:")
        self.runJob("relion_refine --version", "", numberOfMpi=1)
        self.info("Detected version from config: %s"
                  % relion.Plugin.getActiveVersion())

    def _convertCoordinates(self, mic, coordList):
        relion.convert.writeMicCoordinates(
            mic, coordList, self._getMicPos(mic), getPosFunc=self._getPos)

    def _extractMicrograph(self, mic, params):
        """ Extract particles from one micrograph, ignore if the .star
        with the coordinates is not present. """
        self._extractMicrographList([mic], params)

    def _extractMicrographList(self, micList, params):
        workingDir = self.getWorkingDir()
        micsStar = self._getMicsStar(micList)
        og = OpticsGroups.fromImages(self.getInputMicrographs())
        starWriter = relion.convert.createWriter(rootDir=workingDir,
                                                 outputDir=self._getTmpPath(),
                                                 optics=og)
        starWriter.writeSetOfMicrographs(micList, self._getPath(micsStar))

        # convert.writeSetOfMicrographs(micList, self._getPath(micsStar),
        #                               alignType=emcts.ALIGN_NONE,
        #                               preprocessImageRow=self._preprocessMicrographRow)

        partsStar = self._getMicParticlesStar(micList)

        for mic in micList:
            self._convertCoordinates(mic, self.coordDict[mic.getObjId()])
            self._micCoordStarDict[mic.getObjId()] = partsStar

        args = ' --i %s --part_star %s %s' % (micsStar, partsStar, params)

        self.runJob(self._getProgram('relion_preprocess'), args, cwd=workingDir)

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        if self.doNormalize and self.backDiameter > self.boxSize:
            errors.append("Background diameter for normalization should "
                          "be equal or less than the box size.")

        if self.doRescale and self.rescaledSize.get() % 2 == 1:
            errors.append("Only re-scaling to even-sized images is allowed "
                          "in RELION.")

        if self.hasAttribute('saveFloat16') and self.saveFloat16:
            errors.append("MRC float16 format is not yet supported by XMIPP, "
                          "so you cannot use this option.")

        return errors

    def _citations(self):
        return ['Scheres2012b']

    def _summary(self):
        summary = list()
        summary.append("Micrographs source: %s"
                       % self.getEnumText("downsampleType"))
        summary.append("Particle box size: %d" % self.boxSize)

        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.")
        else:
            summary.append("Particles extracted: %d" %
                           self.outputParticles.getSize())

        return summary

    def _methods(self):
        methodsMsgs = []

        if self.getStatus() == STATUS_FINISHED:
            msg = ("A total of %d particles of size %d were extracted"
                   % (self.getOutput().getSize(), self.boxSize))

            if self._micsOther():
                msg += (" from another set of micrographs: %s"
                        % self.getObjectTag('inputMicrographs'))

            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)

            if self.doInvert:
                methodsMsgs.append("Inverted contrast on images.")
            if self._doDownsample():
                methodsMsgs.append("Particles downsampled by a factor of %0.2f."
                                   % self._getDownFactor())

        return methodsMsgs

    # -------------------------- UTILS functions ------------------------------
    def _setupBasicProperties(self):
        # Set sampling rate (before and after doDownsample) and inputMics
        # according to micsSource type
        inputCoords = self.getCoords()
        self.samplingInput = inputCoords.getMicrographs().getSamplingRate()
        self.samplingMics = self.getInputMicrographs().getSamplingRate()
        self.samplingFactor = float(self.samplingMics / self.samplingInput)

        scale = self.getScaleFactor()
        self.debug("Scale: %f" % scale)
        if self.notOne(scale):
            # If we need to scale the box, then we need to scale the coordinates
            getPos = lambda coord: (int(coord.getX() * scale),
                                    int(coord.getY() * scale))
        else:
            getPos = lambda coord: coord.getPosition()
        # Store the function to be used for scaling coordinates
        self._getPos = getPos

    def _getExtractArgs(self):
        # The following parameters are executing 'relion_preprocess' to
        # extract the particles of a given micrographs
        # The following is assumed:
        # - relion_preprocess will be executed from the protocol workingDir
        # - the micrographs (or links) and coordinate files will be in 'extra'
        # - coordinate files have the 'coords.star' suffix
        params = ' --coord_dir "."'
        params += ' --coord_suffix .coords.star'
        params += ' --part_dir "." --extract '
        params += ' --extract_size %d' % self.boxSize

        if Plugin.IS_GT31() and self.saveFloat16:
            params += "--float16 "

        if self.backDiameter <= 0:
            diameter = self.boxSize.get() * 0.75
        else:
            diameter = self.backDiameter.get()

        params += ' --bg_radius %d' % int(diameter/(2 * self._getDownFactor()))

        if self.doInvert:
            params += ' --invert_contrast'

        if self.doNormalize:
            params += ' --norm'

        if self._doDownsample():
            params += ' --scale %d' % self.rescaledSize

        if self.stddevWhiteDust > 0:
            params += ' --white_dust %0.3f' % self.stddevWhiteDust

        if self.stddevBlackDust > 0:
            params += ' --black_dust %0.3f' % self.stddevBlackDust

        return [params]

    def readPartsFromMics(self, micList, outputParts):
        """ Read the particles extract for the given list of micrographs
        and update the outputParts set with new items.
        """
        relionToLocation = relion.convert.relionToLocation
        p = Particle()
        p._rlnOpticsGroup = Integer()
        acq = self.getInputMicrographs().getAcquisition()
        # JMRT: Ideally I would like to disable the whole Acquisition for each
        #       particle row, but the SetOfImages will set it again.
        #       Another option could be to disable in the set, but then in
        #       streaming, other protocols might get the wrong optics info
        pAcq = Acquisition(magnification=acq.getMagnification(),
                           voltage=acq.getVoltage(),
                           amplitudeContrast=acq.getAmplitudeContrast(),
                           sphericalAberration=acq.getSphericalAberration())
        p.setAcquisition(pAcq)

        tmp = self._getTmpPath()
        extra = self._getExtraPath()

        for mic in micList:
            posSet = set()
            coordDict = {self._getPos(c): c
                         for c in self.coordDict[mic.getObjId()]}
            del self.coordDict[mic.getObjId()]

            ogNumber = mic.getAttributeValue('_rlnOpticsGroup', 1)

            partsStar = self.__getMicFile(mic, '_extract.star', folder=tmp)
            partsTable = relion.convert.Table(fileName=partsStar)
            stackFile = self.__getMicFile(mic, '.mrcs', folder=tmp)
            endStackFile = self.__getMicFile(mic, '.mrcs', folder=extra)
            pwutils.moveFile(stackFile, endStackFile)

            for part in partsTable:
                pos = (int(float(part.rlnCoordinateX)),
                       int(float(part.rlnCoordinateY)))

                if pos in posSet:
                    print("Duplicate coordinate at: %s, IGNORED. " % str(pos))
                    coord = None
                else:
                    coord = coordDict.get(pos, None)

                if coord is not None:
                    # scale the coordinates according to particles dimension.
                    coord.scale(self.getBoxScale())
                    p.copyObjId(coord)
                    idx, fn = relionToLocation(part.rlnImageName)
                    p.setLocation(idx, endStackFile)
                    p.setCoordinate(coord)
                    p.setMicId(mic.getObjId())
                    p.setCTF(mic.getCTF())
                    p._rlnOpticsGroup.set(ogNumber)
                    outputParts.append(p)
                    posSet.add(pos)

    def _updateOutputSet(self, outputName, outputSet,
                         state=Set.STREAM_OPEN):
        """ Redefine this method to update optics info. """

        first = getattr(self, '_firstUpdate', True)

        if first:
            og = OpticsGroups.fromImages(outputSet)
            og.updateAll(rlnImagePixelSize=self._getNewSampling(),
                         rlnImageSize=self.getNewImgSize())
            og.toImages(outputSet)

        ProtExtractParticles._updateOutputSet(self, outputName, outputSet,
                                              state=state)
        self._firstUpdate = False

    def _micsOther(self):
        """ Return True if other micrographs are used for extract. """
        return self.downsampleType == OTHER

    def _getDownFactor(self):
        if self.doRescale:
            return float(self.boxSize.get()) / self.rescaledSize.get()
        return 1

    def _doDownsample(self):
        return self.doRescale and self.rescaledSize != self.boxSize

    def notOne(self, value):
        return abs(value - 1) > 0.0001

    def _getNewSampling(self):
        newSampling = self.getInputMicrographs().getSamplingRate()

        if self._doDownsample():
            # Set new sampling, it should be the input sampling of the used
            # micrographs multiplied by the downFactor
            newSampling *= self._getDownFactor()

        return newSampling

    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates or
        Other micrographs. """
        if not self._micsOther():
            return self.inputCoordinates.get().getMicrographs()
        else:
            return self.inputMicrographs.get()

    def getCoords(self):
        return self.inputCoordinates.get()

    def getOutput(self):
        if (self.hasAttribute('outputParticles') and
                self.outputParticles.hasValue()):
            return self.outputParticles
        else:
            return None

    def getCoordSampling(self):
        return

    def getScaleFactor(self):
        """ Returns the scaling factor that needs to be applied to the input
        coordinates to adapt for the input micrographs.
        """
        coordsSampling = self.getCoords().getMicrographs().getSamplingRate()
        micsSampling = self.getInputMicrographs().getSamplingRate()
        return coordsSampling / micsSampling

    def getBoxScale(self):
        """ Computing the sampling factor between input and output.
        We should take into account the differences in sampling rate between
        micrographs used for picking and the ones used for extraction.
        The downsampling factor could also affect the resulting scale.
        """
        f = self.getScaleFactor()
        return f / self._getDownFactor() if self._doDownsample() else f

    def getNewImgSize(self):
        return int(self.rescaledSize if self._doDownsample() else self.boxSize)

    def _getOutputImgMd(self):
        return self._getPath('images.xmd')

    def createParticles(self, item, row):
        particle = relion.convert.rowToParticle(row, readCtf=self._useCTF())
        coord = particle.getCoordinate()
        item.setY(coord.getY())
        item.setX(coord.getX())
        particle.setCoordinate(item)

        item._appendItem = False

    def __getMicFile(self, mic, ext, folder=None):
        """ Return a filename based on the micrograph.
        The filename will be located in the extra folder and with
        the given extension.
        """
        fileName = 'mic_%06d%s' % (mic.getObjId(), ext)
        dirName = folder or self._getExtraPath()
        return os.path.join(dirName, fileName)

    def _useCTF(self):
        return self.ctfRelations.hasValue()

    def _getMicStarFile(self, mic):
        return self.__getMicFile(mic, '.star')

    def _getMicStackFile(self, mic):
        return self.__getMicFile(mic, '.mrcs')

    def _getMicParticlesStar(self, micList):
        """ Return the star files with the particles for this micrographs. """
        return self._getMicsStar(micList).replace('.star', '_particles.star')

    def _getMicPos(self, mic):
        """ Return the corresponding .pos file for a given micrograph. """
        return self.__getMicFile(mic, ".coords.star", folder=self._getTmpPath())

    def _getMicsStar(self, micList):
        return 'micrographs_%05d-%05d.star' % (micList[0].getObjId(),
                                               micList[-1].getObjId())

    def _isStreamOpen(self):
        if self._useCTF():
            ctfStreamOpen = self.ctfRelations.get().isStreamOpen()
        else:
            ctfStreamOpen = False

        return (self.getInputMicrographs().isStreamOpen() or
                ctfStreamOpen or self.getCoords().isStreamOpen())

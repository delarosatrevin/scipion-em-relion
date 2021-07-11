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

from os.path import relpath, abspath

import pyworkflow.protocol.params as params
from pwem.protocols import ProtParticlePickingAuto
from pwem.constants import RELATION_CTF
from pwem.emlib.image import ImageHandler
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD

from relion import Plugin
import relion.convert
from ..constants import *
from .protocol_autopick import ProtRelionAutopickBase


class ProtRelion2Autopick(ProtRelionAutopickBase):
    """ This protocol runs Relion autopicking (version > 3.0).

    This Relion protocol uses the 'relion_autopick' program to pick particles
    from micrographs, either using references (2D averages or 3D volumes)

    The wrapper implementation does not read/write any FOM maps compared to Relion
    """
    _label = 'auto-picking'
    _devStatus = PROD

    # -------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs', important=True,
                      help='Select the input micrographs. '
                           'If using the *Optimize* mode, just a subset of '
                           'micrographs are used to compute the FOM maps. '
                           'If in *Compute* mode, all micrographs will be '
                           'auto-picked.')
        form.addParam('ctfRelations', params.RelationParam,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to the '
                           'input micrographs.')

        # From Relion 3.+, references can be 2D or 3D
        # need to add these parameters
        refCondition = 'referencesType==%s' % REF_AVERAGES
        ref3dCondition = 'referencesType==%s' % REF_VOLUME

        form.addSection('References')

        form.addParam('referencesType', params.EnumParam,
                      choices=['2D', '3D'],
                      default=REF_AVERAGES,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='References',
                      help='The preferred way to autopick is '
                           'by providing 2D references images that were '
                           'obtained by 2D classification. \n'
                           'The Gaussian blob references may be useful to '
                           'kickstart a new data set.')

        # In Relion 3 it is also possible to pass a volume as reference for
        # autopicking

        form.addParam('inputReferences', params.PointerParam,
                      pointerClass='SetOfAverages',
                      condition=refCondition,
                      label='Input references', important=True,
                      help='Input references (SetOfAverages) for auto-pick. \n\n'
                           'Note that the absolute greyscale needs to be correct, \n'
                           'so only use images with proper normalization.')

        form.addParam('inputReferences3D', params.PointerParam,
                      pointerClass='Volume',
                      condition=ref3dCondition,
                      label='Input references', important=True,
                      help='Input volume from which 2D references will be '
                           'made by projection. Note that the absolute '
                           'greyscale needs to be correct, so only use '
                           'maps created by RELION itself from this data set.')

        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      condition=ref3dCondition,
                      label='Symmetry',
                      help="Symmetry point group of the 3D reference. "
                           "Only projections in the asymmetric part of the "
                           "sphere will be generated.")

        form.addParam('angularSamplingDeg', params.EnumParam, default=0,
                      choices=ANGULAR_SAMPLING_LIST,
                      condition=ref3dCondition,
                      label='3D angular sampling (deg)',
                      help="There are only a few discrete angular samplings "
                           "possible because we use the HealPix library to "
                           "generate the sampling of the first two Euler "
                           "angles on the sphere. The samplings are approximate "
                           "numbers and vary slightly over the sphere.\n"
                           "For autopicking, 30 degrees is usually fine enough, "
                           "but for highly symmetrical objects one may need to "
                           "go finer to adequately sample the asymmetric part of "
                           "the sphere.")

        if not Plugin.IS_GT31():
            form.addParam('particleDiameter', params.IntParam, default=-1,
                          label='Mask diameter (A)',
                          help='Diameter of the circular mask that will be applied '
                               'around the templates in Angstroms. When set to a '
                               'negative value, this value is estimated '
                               'automatically from the templates themselves.')

        form.addParam('lowpassFilterRefs', params.IntParam, default=20,
                      label='Lowpass filter references (A)',
                      help='Lowpass filter that will be applied to the '
                           'references before template matching. \n'
                           'Do NOT use very high-resolution templates to '
                           'search your micrographs. \n'
                           'The signal will be too weak at high resolution '
                           'anyway, and you may find Einstein from noise...')

        form.addParam('highpassFilterMics', params.IntParam, default=-1,
                      label='Highpass filter (A)',
                      help='Highpass filter that will be applied to the '
                           'micrographs. This may be useful to get rid of '
                           'background ramps due to uneven ice distributions. '
                           'Give a negative value to skip the highpass '
                           'filter.  Useful values are often in the range '
                           'of 200-400 Angstroms.')

        form.addParam('angularSampling', params.IntParam, default=5,
                      label='Angular sampling (deg)',
                      help='Angular sampling in degrees for exhaustive searches '
                           'of the in-plane rotations for all references.')

        form.addParam('refsHaveInvertedContrast', params.BooleanParam,
                      default=True,
                      label='References have inverted contrast?',
                      help='Set to Yes to indicate that the reference have '
                           'inverted contrast with respect to the particles '
                           'in the micrographs.')

        form.addParam('refsCtfCorrected', params.BooleanParam, default=True,
                      label='Are References CTF corrected?',
                      help='Set to Yes if the references were created with '
                           'CTF-correction inside RELION.\n'
                           'If set to Yes, the input micrographs should contain '
                           'the CTF information.')

        form.addParam('ignoreCTFUntilFirstPeak', params.BooleanParam,
                      default=False,
                      label='Ignore CTFs until first peak?',
                      help='Set this to Yes, only if this option was also used '
                           'to generate the references.')

        form.addSection('Autopicking')

        group = form.addGroup('Autopick')
        group.addParam('pickingThreshold', params.FloatParam, default=0.25,
                       label='Picking threshold:',
                       help='Use lower thresholds to pick more particles '
                            '(and more junk probably)')

        group.addParam('interParticleDistance', params.IntParam, default=-1,
                       label='Minimum inter-particle distance (A):',
                       help='Particles closer together than this distance \n'
                            'will be consider to be a single cluster. \n'
                            'From each cluster, only one particle will be '
                            'picked.')

        group.addParam('maxStddevNoise', params.FloatParam, default=1.1,
                       label='Maximum stddev noise:',
                       help='This is useful to prevent picking in carbon areas, '
                            'or areas with big contamination features. Peaks in '
                            'areas where the background standard deviation in '
                            'the normalized micrographs is higher than this '
                            'value will be ignored. Useful values are probably '
                            'in the range 1.0 to 1.2. Set to -1 to switch off '
                            'the feature to eliminate peaks due to high '
                            'background standard deviations.')

        group.addParam('minAvgNoise', params.FloatParam, default=-999,
                       label='Minimum avg noise:',
                       help='This is useful to prevent picking in carbon areas,'
                            ' or areas with big contamination features. Peaks '
                            'in areas where the background standard deviation '
                            'in the normalized micrographs is higher than this'
                            ' value will be ignored. Useful values are '
                            'probably in the range -0.5 to 0. Set to -999 to '
                            'switch off the feature to eliminate peaks due to '
                            'low average background densities.')

        group = form.addGroup('Computing')
        group.addParam('shrinkFactor', params.FloatParam, default=0,
                       validators=[params.Range(0, 1, "value should be "
                                                      "between 0 and 1. ")],
                       label='Shrink factor',
                       help='This is useful to speed up the calculations, '
                            'and to make them less memory-intensive. The '
                            'micrographs will be downscaled (shrunk) to '
                            'calculate the cross-correlations, and peak '
                            'searching will be done in the downscaled FOM '
                            'maps. When set to 0, the micrographs will de '
                            'downscaled to the lowpass filter of the '
                            'references, a value between 0 and 1 will '
                            'downscale the micrographs by that factor. '
                            'Note that the results will not be exactly '
                            'the same when you shrink micrographs!')

        group.addParam('doGpu', params.BooleanParam, default=True,
                       label='Use GPU acceleration?',
                       help='If set to Yes, the job will try to use GPU '
                            'acceleration.')

        group.addParam('gpusToUse', params.StringParam, default='',
                       label='Which GPUs to use:', condition='doGpu',
                       help='This argument is not necessary. If left empty, '
                            'the job itself will try to allocate available GPU '
                            'resources. You can override the default '
                            'allocation by providing a list of which GPUs '
                            '(0,1,2,3, etc) to use. MPI-processes are '
                            'separated by ":", threads by ",". '
                            'For example: "0,0:1,1:0,0:1,1"')

        form.addParam('extraParams', params.StringParam, default='',
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "relion_autopick" will print a list '
                           'of possible options.')

        self._defineStreamingParams(form)

        form.addParallelSection(threads=0, mpi=4)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.inputStreaming = self.getInputMicrographs().isStreamOpen()

        if self.streamingBatchSize > 0 or self.inputStreaming:
            # If the input is in streaming, follow the base class policy
            # about inserting new steps and discovery new input/output
            self.createOutputStep = self._doNothing
            ProtParticlePickingAuto._insertAllSteps(self)
        else:
            # If not in streaming, then we will just insert a single step to
            # pick all micrographs at once since it is much faster
            self.micDict = {}
            self.micDict, _ = self._loadInputList()

            self._insertFunctionStep('convertInputStep',
                                     self.getInputMicrographs().strId(),
                                     self.getInputReferences().strId())
            nameList = [mic.getMicName() for mic in self.getMicrographList()]
            self._insertFunctionStep('pickMicrographListStep', nameList,
                                     *self._getPickArgs())
            self._insertFunctionStep('createOutputStep')

            # Disable streaming functions:
            self._insertFinalSteps = self._doNothing
            self._stepsCheck = self._doNothing

    def _insertInitialSteps(self):
        # Convert the input micrographs and references to
        # the required Relion star files
        inputRefs = self.getInputReferences()
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.getInputMicrographs().strId(),
                                             inputRefs.strId())
        return [convertId]

    def _doNothing(self, *args):
        pass

    def _loadInputList(self):
        """ This function is re-implemented in this protocol, because it have
         a SetOfCTF as input, so for streaming, we only want to report those
         micrographs for which the CTF is ready.
        """
        micDict, micClose = self._loadMics(self.getInputMicrographs())
        ctfDict, ctfClosed = self._loadCTFs(self.ctfRelations.get())

        # Keep the micrographs that have CTF
        # and set the CTF property for those who have it
        readyMics = dict()

        for micKey, mic in micDict.items():
            if micKey in ctfDict:

                mic.setCTF(ctfDict[micKey])
                readyMics[micKey] = mic

        # Return the updated micDict and the closed status
        return readyMics, micClose and ctfClosed

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self, micsId, refsId):
        pwutils.makePath(self._getExtraPath('DONE'))  # Required to report finished

        inputRefs = self.getInputReferences()
        if self.useInputReferences():
            relion.convert.writeReferences(
                inputRefs, self._getPath('reference_2d'), useBasename=True)
        else:
            ImageHandler().convert(inputRefs, self._getPath('reference_3d.mrc'))

    def getAutopickParams(self):
        # Return the autopicking parameters except for the interactive ones:
        # - threshold
        # - minDistance
        # - maxStd
        params = ' --pickname autopick'
        params += ' --odir "./"'
        params += ' --angpix %0.5f' % self.getInputMicrographs().getSamplingRate()
        params += ' --shrink %0.3f' % self.shrinkFactor

        if not Plugin.IS_GT31():
            params += ' --particle_diameter %d' % self.particleDiameter

        if self.doGpu:
            params += ' --gpu "%s"' % self.gpusToUse

        if self.useInputReferences():
            params += ' --ref %s' % abspath(self._getPath('reference_2d.stk'))
        else:  # 3D reference
            params += ' --ref %s' % abspath(self._getPath('reference_3d.mrc'))
            params += ' --sym %s' % self.symmetryGroup
            params += ' --healpix_order %d' % (int(self.angularSamplingDeg.get()) + 1)

        ps = self.getInputReferences().getSamplingRate()
        params += ' --angpix_ref %0.5f' % ps

        if self.refsHaveInvertedContrast:
            params += ' --invert'

        if self.refsCtfCorrected:
            params += ' --ctf'

        if self.ignoreCTFUntilFirstPeak:
            params += ' --ctf_intact_first_peak'

        params += ' --ang %d' % self.angularSampling
        # Negative values for filters means no-filter
        if self.lowpassFilterRefs > 0:
            params += ' --lowpass %d' % self.lowpassFilterRefs
        if self.highpassFilterMics > 0:
            params += ' --highpass %d' % self.highpassFilterMics

        # Add extra params is any
        params += ' %s' % self.extraParams

        return params

    def _getPickArgs(self):
        basicArgs = self.getAutopickParams()
        threshold = self.pickingThreshold.get()
        interDist = self.interParticleDistance.get()
        maxStd = self.maxStddevNoise.get()
        minAvg = self.minAvgNoise.get()
        return [basicArgs, threshold, interDist, maxStd, minAvg]

    def _pickMicrographsFromStar(self, micStarFile, cwd, params,
                                 threshold, minDistance,
                                 maxStddevNoise, minAvgNoise):
        """ Launch the 'relion_autopick' for micrographs in the inputStarFile.
         If the input set of complete, the star file will contain all the
         micrographs. If working in streaming, it will be only one micrograph.
        """
        params += ' --i %s' % relpath(micStarFile, cwd)
        params += ' --threshold %0.3f' % threshold
        params += ' --min_distance %0.3f' % minDistance
        params += ' --max_stddev_noise %0.3f' % maxStddevNoise
        params += ' --min_avg_noise %0.3f' % minAvgNoise

        program = self._getProgram('relion_autopick')

        self.runJob(program, params, cwd=cwd)

    # -------------------------- STEPS functions -------------------------------
    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        outputCoordinatesName = 'outputCoordinates'
        coordSet = self._createSetOfCoordinates(micSet)
        self.readCoordsFromMics(None, micSet, coordSet)

        self._defineOutputs(**{outputCoordinatesName: coordSet})
        self._defineSourceRelation(self.getInputMicrographsPointer(),
                                   coordSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        if self.useInputReferences():
            if not Plugin.IS_GT31() and (self.particleDiameter > self.getInputDimA()):
                errors.append('Particle diameter (%d) can not be greater than '
                              'size (%d)' % (self.particleDiameter,
                                             self.getInputDimA()))
            if self.getInputReferences().isOddX():
                errors.append("Relion only works with even values for the "
                              "average dimensions!")
        else:
            if self.particleDiameter <= 0:
                errors.append('When using Gaussian blobs, you need to specify '
                              'the particles diameter manually. ')

        if self.ctfRelations.get() is None and self.refsCtfCorrected:
            errors.append("References CTF corrected parameter must be set to "
                          "False or set ctf relations.")

        return errors

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Scheres2015']

    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append("%s: User picked %d particles with a particle "
                               "size of %d px."
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize()))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    # -------------------------- UTILS functions -------------------------------
    def useInputReferences(self):
        return self.referencesType == REF_AVERAGES

    def getInputDimA(self):
        """ Return the dimension of input references in A. """
        inputRefs = self.getInputReferences()
        if inputRefs is None:
            return None
        else:
            return inputRefs.getXDim() * inputRefs.getSamplingRate()

    def getBoxSize(self):
        """ Return a reasonable box-size in pixels. """
        inputRefs = self.getInputReferences()
        inputMics = self.getInputMicrographs()
        micsSampling = inputMics.getSamplingRate()

        if inputRefs is None:
            if Plugin.IS_GT31():
                boxSize = 128
            else:
                boxSize = int(self.particleDiameter.get() * 1.25 / micsSampling)
        else:
            # Scale boxsize if the pixel size of the references is not the same
            # of the micrographs
            scale = inputRefs.getSamplingRate() / micsSampling
            boxSize = int(inputRefs.getXDim() * scale)

        if boxSize % 2 == 1:
            boxSize += 1  # Use even box size for relion

        return boxSize

    def getInputReferences(self):
        if self.useInputReferences():
            return self.inputReferences.get()
        else:
            return self.inputReferences3D.get()

    def getMicrographList(self):
        """ Return the list of micrographs (either a subset or the full set)
        that will be used for optimizing the parameters or the picking.
        """
        # Use all micrographs only when going for the full picking
        inputMics = self.getInputMicrographs()

        return inputMics

    def __getMicListPrefix(self, micList):
        n = len(micList)
        if n == 0:
            raise Exception("Empty micrographs list!")
        micsPrefix = 'mic_%06d' % micList[0].getObjId()
        if n > 1:
            micsPrefix += "-%06d" % micList[-1].getObjId()
        return micsPrefix

    def _getMicStarFile(self, micList):
        return self._getTmpPath(self.__getMicListPrefix(micList) + '.star')

    def _createTmpMicsDir(self, micList):
        """ Create a temporary path to work with a list of micrographs. """
        micsDir = self._getTmpPath(self.__getMicListPrefix(micList))
        pwutils.makePath(micsDir)
        return micsDir

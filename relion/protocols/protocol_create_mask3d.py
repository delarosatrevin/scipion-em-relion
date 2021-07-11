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

import math

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtCreateMask3D
from pwem.objects import VolumeMask

import relion.convert as convert
from ..constants import MASK_AND


class ProtRelionCreateMask3D(ProtCreateMask3D):
    """ This protocols creates a 3D mask using Relion.
    The mask is created from a 3d volume or by comparing two input volumes.
    """
    _label = 'create 3d mask'
    _devStatus = PROD

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Mask generation')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass="Volume",
                      label="Input volume",
                      help="Select the volume that will be used to create the mask")
        form.addParam('initialLowPassFilterA', params.FloatParam,
                      default=15.,
                      label='Lowpass filter map by (A)',
                      help='Lowpass filter that will be applied to the input map, '
                           'prior to binarization. To calculate solvent masks, a '
                           'lowpass filter of 15-20A may work well.')
        # TODO: add wizard
        form.addParam('threshold', params.FloatParam, default=0.02,
                      label='Initial binarisation threshold',
                      help="This threshold is used to make an initial binary "
                           "mask from the average of the two unfiltered "
                           "half-reconstructions. If you don't know what "
                           "value to use, display one of the unfiltered "
                           "half-maps in a 3D surface rendering viewer and "
                           "find the lowest threshold that gives no noise "
                           "peaks outside the reconstruction.")

        form.addParam('doCompare', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Compare with another volume to produce a mask?',
                      help='Logical comparison of two input volumes to produce a mask')

        form.addParam('inputVolume2', params.PointerParam,
                      pointerClass="Volume",
                      condition='doCompare',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Input volume (second)",
                      help="Select the volume that will be compared to the first one")

        form.addParam('operation', params.EnumParam, default=MASK_AND,
                      condition='doCompare',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Operation',
                      choices=['AND', 'OR', 'AND_NOT', 'OR_NOT'],
                      help='*AND*: Pixels in the initial mask will be one if the input '
                           'AND the second volume are above the threshold value.\n'
                           '*OR*: Pixels in the initial mask will be one if the input '
                           'OR the second volume are above the threshold value.\n'
                           '*AND_NOT*: pixels in the initial mask will be one if '
                           'the input is above the threshold AND the second volume '
                           'is below it.\n*OR_NOT*: pixels in the initial mask will be '
                           'one if the input is above the threshold OR the second '
                           'volume is below it.')

        form.addParam('extend', params.IntParam, default=3,
                      label='Extend binary mask by (px)',
                      help='The initial binary mask is extended this number of '
                           'pixels in all directions.')

        form.addParam('edge', params.IntParam, default=3,
                      label='Add a soft-edge (px)',
                      help='The extended binary mask is further extended with '
                           'a raised-cosine soft edge of the specified width.')

        form.addParam('doInvert', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Invert final mask',
                      help='Invert the final mask')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.maskFile = self._getExtraPath('mask.mrc')
        self._insertFunctionStep('convertInputStep', self.inputVolume.get().getObjId())
        self._insertFunctionStep('createMaskStep')
        self._insertFunctionStep('createOutputStep')
    
    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, volId):
        self.inputVolFn = convert.convertBinaryVol(
            self.inputVolume.get(), self._getTmpPath())

        if self.doCompare:
            self.inputVol2Fn = convert.convertBinaryVol(
                self.inputVolume2.get(), self._getTmpPath())

    def createMaskStep(self):
        argsDict = {'--i ': self.inputVolFn,
                    '--ini_threshold ': self.threshold.get(),
                    '--extend_inimask ': self.extend.get(),
                    '--width_soft_edge ': self.edge.get(),
                    '--angpix ': self.inputVolume.get().getSamplingRate()
                    }
        if self.initialLowPassFilterA.get() != -1:
            argsDict['--lowpass '] = self.initialLowPassFilterA.get()

        args = ' --o %s ' % self.maskFile
        args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])

        if self.doCompare:
            op = self.operation.get()
            opStr = ['--and', '--or', '--and_not', '--or_not'][op]
            args += ' %s %s' % (opStr, self.inputVol2Fn)

        if self.doInvert:
            args += ' --invert'

        if self.numberOfThreads > 1:
            args += ' --j %d' % self.numberOfThreads

        self.runJob("relion_mask_create", args)

        return [self.maskFile]

    def createOutputStep(self):
        volMask = VolumeMask()
        volMask.setFileName(self.maskFile)
        volMask.setSamplingRate(self.inputVolume.get().getSamplingRate())

        self._defineOutputs(outputMask=volMask)
        self._defineSourceRelation(self.inputVolume, self.outputMask)
        
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.doCompare:
            pixVol1 = self.inputVolume.get().getSamplingRate()
            pixVol2 = self.inputVolume2.get().getSamplingRate()

            if not math.isclose(pixVol1, pixVol2, abs_tol=0.001):
                errors.append('Pixel size should be the same for both input volumes.')

        return errors

    def _summary(self):
        messages = [
            "Created a mask from input volume using threshold %0.3f"
            % self.threshold.get(),
            "*Mask processing*",
            "   Extend by %d pixels" % self.extend,
            "   Apply soft edge of %d pixels" % self.edge]
        if self.doInvert:
            messages.append("   Inverted")

        return messages

    def _citations(self):
        return []

    def _methods(self):
        messages = [
            "*Mask creation*",
            "We processed the volume %s." % self.getObjectTag('inputVolume'),
            "We binarized it at threshold of %0.3f. " % self.threshold,
            "We extended binary mask by %d voxels." % self.extend,
            "And, we smoothed it by applying a soft edge of %d voxels."
            % self.edge.get()
        ]
        if self.doInvert:
            messages.append("We inverted the mask. ")
            messages.append("And, we smoothed it by applying a soft edge of "
                            "%d voxels." % self.edge)
        if self.hasAttribute('outputMask'):
            messages.append('We refer to the output mask as %s.'
                            % self.getObjectTag('outputMask'))

        return messages

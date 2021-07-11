# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

from pyworkflow.object import String, Integer
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, BooleanParam,
                                        IntParam, LabelParam, LEVEL_ADVANCED)
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtOperateParticles

import relion.convert as convert
from relion import Plugin
from .protocol_base import ProtRelionBase


class ProtRelionSubtract(ProtOperateParticles, ProtRelionBase):
    """ Signal subtraction protocol of Relion.

    Subtract volume projections from the experimental particles.
    The particles must have projection alignment in order to
    properly generate volume projections.
    """
    _label = 'subtract projection'
    _devStatus = PROD

    def _initialize(self):
        self._createFilenameTemplates()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'input_star': self._getExtraPath('input_particles.star'),
                  'output_star': self._getExtraPath('particles_subtracted.star')
                  }
        self._updateFilenamesDict(myDict)
    
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('relionInput', BooleanParam,
                      default=True, important=True,
                      label="Start from Relion protocol?",
                      help="Set to Yes if you wish to use as input "
                           "a Relion protocol. Otherwise set it to No")
        form.addParam('inputProtocol', PointerParam,
                      important=True,
                      pointerClass='ProtRelionRefine3D, ProtRelionClassify3D',
                      label="Input Relion protocol",
                      condition="relionInput",
                      help="Select the 3D refinement/classification run which "
                           "you want to use for subtraction. It will use the "
                           "maps from this run for the subtraction.")
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input map to be projected",
                      important=True,
                      condition="not relionInput",
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
        form.addParam('useAll', BooleanParam, default=True,
                      label="Use all particles from input protocol?",
                      condition="relionInput",
                      help="If No, then you need to provide a subset of "
                           "particles below.")

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not useAll',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles subset",
                      help='Select the particles which are a SUBSET of the '
                           'input protocol provided above.')

        form.addParam('inputParticlesAll', PointerParam,
                      pointerClass='SetOfParticles',
                      condition="not relionInput",
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the input particles.')

        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask of the signal to keep',
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is black (0) and the density you "
                           "wish to keep is white (1).\n"
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to KEEP.*")

        if Plugin.IS_GT31():
            form.addParam('saveFloat16', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Write output in float16?",
                          help="Relion can write output images in float16 "
                               "MRC (mode 12) format to save disk space. "
                               "By default, float32 format is used.")

        form.addSection('Centering')
        form.addParam('help1', LabelParam,
                      condition="not relionInput",
                      label="This section is only used if "
                            "starting from Relion input.")
        form.addParam('centerOnMask', BooleanParam, default=True,
                      label="Do center subtracted images on mask?",
                      condition="relionInput",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the "
                           "center-of-mass of the input mask.")
        form.addParam('centerOnCoord', BooleanParam, default=False,
                      condition='(not centerOnMask) and relionInput',
                      label="Do center on my coordinates?",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the x,y,z "
                           "coordinates below. The unit is pixel, not "
                           "angstrom. The origin is at the center of the box, "
                           "not at the corner.")

        line = form.addLine('Center coordinate (px)',
                            condition='centerOnCoord and relionInput',
                            help='Coordinate of the 3D center (in pixels).')
        line.addParam('cX', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='X')
        line.addParam('cY', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='Y')
        line.addParam('cZ', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='Z')

        form.addParam('newBoxSize', IntParam, default=-1,
                      condition="relionInput",
                      label="New box size",
                      help="Provide a non-negative value to re-window the "
                           "subtracted particles in a smaller box size.")

        form.addSection(label='CTF')
        form.addParam('help', LabelParam,
                      condition="relionInput",
                      label="This section is only used if "
                            "NOT starting from Relion protocol.")
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?',
                      condition="not relionInput",
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      label='Ignore CTFs until first peak?',
                      condition="not relionInput",
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')

        form.addParallelSection(threads=0, mpi=1)
    
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.isRelionInput = self.relionInput.get()
        self._initialize()

        if not self.useAll or not self.isRelionInput:
            self._insertFunctionStep('convertInputStep')

        self._insertFunctionStep('subtractStep')
        self._insertFunctionStep('createOutputStep')
    
    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ Write the input images as a Relion star file. """
        if self.isRelionInput:
            imgSet = self.inputParticles.get()
        else:
            imgSet = self.inputParticlesAll.get()

        convert.writeSetOfParticles(
            imgSet, self._getFileName('input_star'),
            outputDir=self._getExtraPath(), alignType=ALIGN_PROJ)

    def subtractStep(self):
        if self.isRelionInput:
            self.subtractStepRelion()
        else:
            self.subtractStepNoRelion()

    def subtractStepNoRelion(self):
        volume = self.inputVolume.get()
        volFn = convert.convertBinaryVol(volume,
                                         self._getExtraPath())
        params = ' --i %s --subtract_exp' % volFn
        params += ' --angpix %0.3f' % volume.getSamplingRate()
        params += self._convertMask(resize=False, invert=True)


        if self.doCTF:
            params += ' --ctf'
            if self.ignoreCTFUntilFirstPeak:
                params += ' --ctf_intact_first_peak'
            if self._getInputParticles().isPhaseFlipped():
                params += ' --ctf_phase_flip'

        params += ' --ang %s  --o %s' % (
            self._getFileName('input_star'),
            self._getFileName('output_star').replace(".star", ""))

        self.runJob('relion_project', params)

    def subtractStepRelion(self):
        inputProt = self.inputProtocol.get()
        inputProt._initialize()
        fnOptimiser = inputProt._getFileName('optimiser',
                                             iter=inputProt._lastIter())
        params = " --i %s --o %s --new_box %s" % (fnOptimiser,
                                                  self._getExtraPath(),
                                                  self.newBoxSize.get())

        if not self.useAll:
            params += " --data %s" % self._getFileName('input_star')
        if self.centerOnMask:
            params += " --recenter_on_mask"
        elif self.centerOnCoord:
            params += " --center_x %d --center_y %d --center_z %d" % (
                self.cX, self.cY, self.cZ)

        if Plugin.IS_GT31() and self.saveFloat16:
            params += " --float16"

        params += self._convertMask()
        self.runJob(self._getProgram('relion_particle_subtract'), params)

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        outImgSet = self._createSetOfParticles()
        outImgsFn = self._getFileName('output_star')
        outImgSet.copyInfo(imgSet)
        outImgSet.setAlignmentProj()

        px = imgSet.getSamplingRate()
        self.reader = convert.createReader(alignType=ALIGN_PROJ,
                                           pixelSize=px)
        mdIter = convert.Table.iterRows('particles@' + outImgsFn)
        outImgSet.copyItems(imgSet, doClone=False,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=mdIter)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if not self.useAll:
            self._validateDim(self.inputParticles(),
                              self._getInputParticles().getXDim(),
                              errors, 'Input particles subset',
                              'Input particles from 3D protocol')
        if self.numberOfMpi > 1 and (not self.relionInput.get()):
            errors.append("Use of several CPUs when input is not relion "
                          "protocol is not supported")
        if self.hasAttribute('saveFloat16') and self.saveFloat16:
            errors.append("MRC float16 format is not yet supported by XMIPP, "
                          "so you cannot use this option.")
        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output is not ready yet.")
        else:
            summary.append('Projections of the masked input volume were '
                           'subtracted from original particles.')

        return summary
    
    # -------------------------- UTILS functions ------------------------------
    def _updateItem(self, particle, row):
        if self.isRelionInput:
            # FIXME: check if other attrs need saving
            particle._rlnRandomSubset = Integer(row.rlnRandomSubset)
            self.reader.setParticleTransform(particle, row)
        particle._rlnImageOriginalName = String(row.rlnImageOriginalName)
        newFn = row.rlnImageName
        newLoc = convert.relionToLocation(newFn)
        particle.setLocation(newLoc)

    def _getInputParticles(self):
        if self.isRelionInput:
            inputProt = self.inputProtocol.get()
            return inputProt.outputParticles
        else:
            return self.inputParticlesAll.get()

    def _convertMask(self, invert=False, resize=True):
        tmp = self._getTmpPath()
        if resize:
            newDim = self._getInputParticles().getXDim()
            newPix = self._getInputParticles().getSamplingRate()
        else:
            newDim = None
            newPix = None
        maskFn = convert.convertMask(self.refMask.get(),
                                     tmp, newPix, newDim, invert=invert)
        return ' --mask %s' % maskFn

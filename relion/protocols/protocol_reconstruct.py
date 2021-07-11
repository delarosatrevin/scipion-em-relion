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

from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam,
                                        EnumParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.constants import PROD
from pwem.objects import Volume
from pwem.protocols import ProtReconstruct3D
from pwem.constants import ALIGN_PROJ

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionReconstruct(ProtReconstruct3D, ProtRelionBase):
    """ This protocol reconstructs a volume using Relion.

    Reconstruct a volume from a given set of particles.
    The alignment parameters will be converted to a Relion star file
    and used as direction projections to reconstruct.
    """
    _label = 'reconstruct'
    _devStatus = PROD

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group",
                      help='See [[https://relion.readthedocs.io/'
                           'en/latest/Reference/Conventions.html#symmetry]'
                           '[Relion Symmetry]] page for a description '
                           'of the symmetry format accepted by Relion')
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).')
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor")
        form.addParam('subset', EnumParam, default=0,
                      choices=['all', 'half1', 'half2'],
                      display=EnumParam.DISPLAY_HLIST,
                      label='Subset to reconstruct',
                      help='Subset of images to consider.')
        form.addParam('classNum', IntParam, default=-1,
                      label='Use only this class',
                      help='Consider only this class (-1: use all classes)')
        
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *relion_reconstruct* program. '
                           'Address to Relion to see full list of options.')
        form.addSection('CTF')
        form.addParam('doCTF', BooleanParam, default=False,
                      label='Apply CTF correction?')
        form.addParam('ctfIntactFirstPeak', BooleanParam, default=False,
                      condition='doCTF',
                      label='Leave CTFs intact until first peak?')
        
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertReconstructStep()
        self._insertFunctionStep('createOutputStep')

    def _insertReconstructStep(self):
        imgSet = self.inputParticles.get()

        params = ' --i %s' % self._getFileName('input_particles')
        params += ' --o %s' % self._getFileName('output_volume')
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angpix %0.5f' % imgSet.getSamplingRate()
        params += ' --maxres %0.3f' % self.maxRes.get()
        params += ' --pad %0.3f' % self.pad.get()

        subset = -1 if self.subset.get() == 0 else self.subset
        params += ' --subset %d' % subset
        params += ' --class %d' % self.classNum.get()

        if self.doCTF:
            params += ' --ctf'
            if self.ctfIntactFirstPeak:
                params += ' --ctf_intact_first_peak'

            if imgSet.isPhaseFlipped():
                params += ' --ctf_phase_flipped'

        if self.extraParams.hasValue():
            params += " " + self.extraParams.get()

        self._insertFunctionStep('reconstructStep', params)

    # -------------------------- STEPS functions ------------------------------
    def reconstructStep(self, params):
        self._runProgram('relion_reconstruct', params)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'output_volume': self._getExtraPath('output_volume.mrc')
            }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_particles')

        # Pass stack file as None to avoid write the images files
        convert.writeSetOfParticles(imgSet, imgStar,
                                    outputDir=self._getTmpPath(),
                                    alignType=ALIGN_PROJ)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume = Volume()
        volume.setFileName(self._getFileName('output_volume'))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputParticles, volume)
    
    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
        else:
            summary.append("Output volume has been reconstructed.")

        return summary

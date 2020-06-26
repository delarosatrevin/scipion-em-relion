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
from emtable import Table

from pyworkflow.object import String, Float
import pwem
from pwem.objects import SetOfClasses2D
from pwem.protocols import ProtClassify2D

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionClassify2D(ProtRelionBase, ProtClassify2D):
    """ This protocol runs Relion 2D classification."""

    _label = '2D classification'
    IS_2D = True
    OUTPUT_TYPE = SetOfClasses2D
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'

    # --------------------------- INSERT steps functions ----------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment:
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * self._getSamplingFactor()
            args['--psi_step'] = self.inplaneAngularSamplingDeg.get() * self._getSamplingFactor()

            if self.IS_GT30() and self.allowCoarserSampling:
                args['--allow_coarser_sampling'] = ''

        else:
            args['--skip_align'] = ''

    # --------------------------- STEPS functions -----------------------------
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses2D from a given iteration. """
        classLoader = convert.ClassesLoader(self, pwem.ALIGN_2D)
        classLoader.fillClassesFromIter(clsSet, iteration)

    def createOutputStep(self):
        partSet = self.inputParticles.get()       
        
        classes2D = self._createSetOfClasses2D(partSet)
        self._fillClassesFromIter(classes2D, self._lastIter())
        
        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(self.inputParticles, classes2D)
        
    # --------------------------- INFO functions ------------------------------
    def _validateNormal(self):
        return []
    
    def _validateContinue(self):
        errors = []
        continueRun = self.continueRun.get()
        continueRun._initialize()
        lastIter = continueRun._lastIter()
        
        if self.continueIter.get() == 'last':
            continueIter = lastIter
        else:
            continueIter = int(self.continueIter.get())
        
        if continueIter > lastIter:
            errors += ["You can continue only from the iteration %01d or less" % lastIter]
        
        return errors
    
    def _summaryNormal(self):
        summary = list()
        summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
        summary.append("Classified into *%d* classes." % self.numberOfClasses)
        summary.append("Output set: %s" % self.getObjectTag('outputClasses'))
        
        return summary
    
    def _summaryContinue(self):
        summary = list()
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        
        return summary
    
    def _methods(self):
        methods = ''
        if hasattr(self, 'outputClasses'):
            methods += "We classified input particles %s (%d items) " % (
                self.getObjectTag('inputParticles'),
                self.inputParticles.get().getSize())
            methods += "into %d classes using Relion Classify2d. " % self.numberOfClasses
            methods += 'Output classes: %s' % self.getObjectTag('outputClasses')
        return [methods]
    

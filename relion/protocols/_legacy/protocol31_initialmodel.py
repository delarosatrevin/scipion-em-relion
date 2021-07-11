# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
# * [2] SciLifeLab, Stockholm University
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

from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtInitialVolume
from pwem.objects import Volume
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, IntParam,
                                        EnumParam, StringParam,
                                        BooleanParam,
                                        LEVEL_ADVANCED)

import relion
import relion.convert as convert
from ..protocol_base import ProtRelionBase


class ProtRelionInitialModel(ProtInitialVolume, ProtRelionBase):
    """ This protocols creates a 3D initial model using Relion.

    Generate a 3D initial model _de novo_ from 2D particles using
    Relion Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = '3D initial model'
    IS_CLASSIFY = False
    IS_3D_INIT = True
    IS_2D = False
    CHANGE_LABELS = ['rlnChangesOptimalOrientations',
                     'rlnChangesOptimalOffsets']

    def __init__(self, **args):
        ProtRelionBase.__init__(self, **args)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()
        if not self.doContinue:
            self.continueRun.set(None)
        self.maskZero = False
        self.copyAlignment = False
        self.hasReferenceCTFCorrected = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineConstants()

        self.IS_3D = not self.IS_2D
        form.addSection(label='Input')
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                           'run of type *%s* class and most of the input '
                           'parameters will be taken from it.'
                           % self.getClassName())
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not doContinue',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('maskDiameterA', IntParam, default=-1,
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a '
                           'soft circular mask with this <diameter>. '
                           'Make sure this diameter is not set too small '
                           'because that may mask away part of the signal! If '
                           'set to a value larger than the image size no '
                           'masking will be performed.\n\n'
                           'The same diameter will also be used for a '
                           'spherical mask of the reference structures if no '
                           'user-provided mask is specified.')
        form.addParam('continueRun', PointerParam,
                      pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', StringParam, default='last',
                      condition='doContinue',
                      label='Continue from iteration',
                      help='Select from which iteration do you want to '
                           'continue. If you use *last*, then the last '
                           'iteration will be used. Otherwise, a valid '
                           'iteration number should be provided.')

        self.addSymmetry(form)

        form.addSection(label='CTF')
        form.addParam('continueMsg', LabelParam, default=True,
                      condition='doContinue',
                      label='CTF parameters are not available in continue mode')
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?', condition='not doContinue',
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('haveDataBeenPhaseFlipped', LabelParam,
                      condition='not doContinue',
                      label='Have data been phase-flipped?      '
                            '(Don\'t answer, see help)',
                      help='The phase-flip status is recorded and managed by '
                           'Scipion. \n In other words, when you import or '
                           'extract particles, \nScipion will record whether '
                           'or not phase flipping has been done.\n\n'
                           'Note that CTF-phase flipping is NOT a necessary '
                           'pre-processing step \nfor MAP-refinement in '
                           'RELION, as this can be done inside the internal\n'
                           'CTF-correction. However, if the phases have been '
                           'flipped, the program will handle it.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      label='Ignore CTFs until first peak?',
                      condition='not doContinue',
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')
        form.addParam('doCtfManualGroups', BooleanParam, default=False,
                      label='Do manual grouping ctfs?',
                      condition='not doContinue',
                      help='Set this to Yes the CTFs will grouping manually.')
        form.addParam('defocusRange', FloatParam, default=1000,
                      label='Defocus range for group creation (in Angstroms)',
                      condition='doCtfManualGroups and not doContinue',
                      help='Particles will be grouped by defocus.'
                           'This parameter is the bin for a histogram.'
                           'All particles assigned to a bin form a group')
        form.addParam('numParticles', FloatParam, default=10,
                      label='minimum size for defocus group',
                      condition='doCtfManualGroups and not doContinue',
                      help='If defocus group is smaller than this value, '
                           'it will be expanded until number of particles '
                           'per defocus group is reached')

        form.addSection('Optimisation')
        form.addParam('numberOfClasses', IntParam, default=1,
                      condition='not doContinue',
                      label='Number of classes',
                      help='The number of classes (K) for a multi-reference '
                           'ab initio SGD refinement. These classes will be '
                           'made in an unsupervised manner, starting from a '
                           'single reference in the initial iterations of '
                           'the SGD, and the references will become '
                           'increasingly dissimilar during the in between '
                           'iterations.')

        form.addParam('doFlattenSolvent', BooleanParam, default=True,
                      condition='not doContinue',
                      label='Flatten and enforce non-negative solvent?',
                      help='If set to Yes, the job will apply a spherical '
                           'mask and enforce all values in the reference '
                           'to be non-negative.')

        form.addParam('symmetryGroup', StringParam, default='c1',
                      condition='not doContinue',
                      label="Symmetry",
                      help='SGD sometimes works better in C1. If you make an '
                           'initial model in C1 but want to run Class3D/Refine3D '
                           'with a higher point group symmetry, the reference model '
                           'must be rotated to conform the symmetry convention. '
                           'You can do this by the relion_align_symmetry command.')

        group = form.addGroup('Sampling')
        group.addParam('angularSamplingDeg', EnumParam, default=1,
                       choices=relion.ANGULAR_SAMPLING_LIST,
                       label='Initial angular sampling (deg)',
                       help='There are only a few discrete angular samplings'
                            ' possible because we use the HealPix library to'
                            ' generate the sampling of the first two Euler '
                            'angles on the sphere. The samplings are '
                            'approximate numbers and vary slightly over '
                            'the sphere.')
        group.addParam('offsetSearchRangePix', FloatParam, default=6,
                       label='Offset search range (pix)',
                       help='Probabilities will be calculated only for '
                            'translations in a circle with this radius (in '
                            'pixels). The center of this circle changes at '
                            'every iteration and is placed at the optimal '
                            'translation for each image in the previous '
                            'iteration.')
        group.addParam('offsetSearchStepPix', FloatParam, default=2,
                       label='Offset search step (pix)',
                       help='Translations will be sampled with this step-size '
                            '(in pixels). Translational sampling is also done '
                            'using the adaptive approach. Therefore, if '
                            'adaptive=1, the translations will first be '
                            'evaluated on a 2x coarser grid.')

        form.addSection(label='SGD')
        self._defineSGD3(form)

        form.addParam('sgdNoiseVar', IntParam, default=-1,
                      condition='not doContinue',
                      expertLevel=LEVEL_ADVANCED,
                      label='Increased noise variance half-life',
                      help='When set to a positive value, the initial '
                           'estimates of the noise variance will internally '
                           'be multiplied by 8, and then be gradually '
                           'reduced, having 50% after this many particles '
                           'have been processed. By default, this option '
                           'is switched off by setting this value to a '
                           'negative number. In some difficult cases, '
                           'switching this option on helps. In such cases, '
                           'values around 1000 have found to be useful. '
                           'Change the factor of eight with the additional '
                           'argument *--sgd_sigma2fudge_ini*')

        form.addSection('Compute')
        self._defineComputeParams(form)

        form.addParam('extraParams', StringParam, default='',
                      label='Additional arguments',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program, e.g: \n"
                           "--dont_combine_weights_via_disc\n"
                           "--verb 1\n"
                           "--pad 2\n")

        form.addParallelSection(threads=1, mpi=3)

    def _defineSGD3(self, form):
        """ Define SGD parameters for Relion version 3. """
        group = form.addGroup('Iterations')
        group.addParam('numberOfIterInitial', IntParam, default=50,
                       label='Number of initial iterations',
                       help='Number of initial SGD iterations, at which the '
                            'initial resolution cutoff and the initial subset '
                            'size will be used, and multiple references are '
                            'kept the same. 50 seems to work well in many '
                            'cases. Increase if the correct solution is not '
                            'found.')
        group.addParam('numberOfIterInBetween', IntParam, default=200,
                       label='Number of in-between iterations',
                       help='Number of SGD iterations between the initial and '
                            'final ones. During these in-between iterations, '
                            'the resolution is linearly increased, together '
                            'with the mini-batch or subset size. In case of a '
                            'multi-class refinement, the different references '
                            'are also increasingly left to become dissimilar. '
                            '200 seems to work well in many cases. Increase '
                            'if multiple references have trouble separating, '
                            'or the correct solution is not found.')
        group.addParam('numberOfIterFinal', IntParam, default=50,
                       label='Number of final iterations',
                       help='Number of final SGD iterations, at which the '
                            'final resolution cutoff and the final subset '
                            'size will be used, and multiple references are '
                            'left dissimilar. 50 seems to work well in many '
                            'cases. Perhaps increase when multiple reference '
                            'have trouble separating.')
        group.addParam('writeIter', IntParam, default=10,
                       expertLevel=LEVEL_ADVANCED,
                       label='Write-out frequency (iter)',
                       help='Every how many iterations do you want to write the '
                            'model to disk. Negative value means only write '
                            'out model after entire iteration.')

        line = form.addLine('Resolution (A)',
                            help='This is the resolution cutoff (in A) that '
                                 'will be applied during the initial and final '
                                 'SGD iterations. 35A and 15A respectively '
                                 'seems to work well in many cases.')
        line.addParam('initialRes', IntParam, default=35, label='Initial')
        line.addParam('finalRes', IntParam, default=15, label='Final')

        line = form.addLine('Mini-batch size',
                            help='The number of particles that will be processed '
                                 'during the initial and final iterations. \n\n'
                                 'For initial, 100 seems to work well in many '
                                 'cases. Lower values may result in wider '
                                 'searches of the energy landscape, but possibly '
                                 'at reduced resolutions. \n\n'
                                 'For final, 300-500 seems to work well in many '
                                 'cases. Higher values may result in increased '
                                 'resolutions, but at increased computational '
                                 'costs.')
        line.addParam('initialBatch', IntParam, default=100, label='Initial')
        line.addParam('finalBatch', IntParam, default=500, label='Final')

    def addSymmetry(self, container):
        pass

    # -------------------------- INSERT steps functions -----------------------

    # -------------------------- STEPS functions ------------------------------
    def _getVolumes(self):
        """ Return the list of volumes generated.
        The number of volumes in the list will be equal to
        the number of classes requested by the user in the protocol. """
        # Provide 1 as default value for making it backward compatible
        k = self.getAttributeValue('numberOfClasses', 1)
        pixelSize = self._getInputParticles().getSamplingRate()
        lastIter = self._lastIter()
        volumes = []

        for i in range(1, k + 1):
            vol = Volume(self._getExtraPath('relion_it%03d_class%03d.mrc')
                         % (lastIter, i))
            vol.setSamplingRate(pixelSize)
            volumes.append(vol)

        return volumes

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        volumes = self._getVolumes()

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet, self._lastIter())

        if len(volumes) > 1:
            output = self._createSetOfVolumes()
            output.setSamplingRate(imgSet.getSamplingRate())
            for vol in volumes:
                output.append(vol)
            self._defineOutputs(outputVolumes=output)
        else:
            output = volumes[0]
            self._defineOutputs(outputVolume=output)

        self._defineSourceRelation(self.inputParticles, output)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

    # -------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        errors = []
        return errors

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
            errors += ["The iteration from you want to continue must be "
                       "%01d or less" % lastIter]

        return errors

    def _summaryNormal(self):
        summary = []
        it = self._lastIter() or -1
        if it >= 1:
            table = Table(fileName=self._getFileName('model', iter=it),
                          tableName='model_general')
            row = table[0]
            resol = float(row.rlnCurrentResolution)
            summary.append("Current resolution: *%0.2f*" % resol)
        return summary

    def _summaryContinue(self):
        summary = ["Continue from iteration %01d" % self._getContinueIter()]
        return summary

    # -------------------------- UTILS functions ------------------------------
    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--o': self._getExtraPath('relion'),
                     '--oversampling': '1'
                     })

        if self.doFlattenSolvent:
            args['--flatten_solvent'] = ''
        if not self.doContinue:
            args.update({'--sym': self.symmetryGroup.get()})
        args['--pad'] = 1 if self.skipPadding else 2
        if self.skipGridding:
            args['--skip_gridding'] = ''

        self._setSGDArgs(args)
        self._setSamplingArgs(args)

    def _setSGDArgs(self, args):
        args['--sgd'] = ''
        args['--sgd_ini_iter'] = self.numberOfIterInitial.get()
        args['--sgd_inbetween_iter'] = self.numberOfIterInBetween.get()
        args['--sgd_fin_iter'] = self.numberOfIterFinal.get()
        args['--sgd_write_iter'] = self.writeIter.get()
        args['--sgd_ini_resol'] = self.initialRes.get()
        args['--sgd_fin_resol'] = self.finalRes.get()
        args['--sgd_ini_subset'] = self.initialBatch.get()
        args['--sgd_fin_subset'] = self.finalBatch.get()
        args['--K'] = self.numberOfClasses.get()

        if not self.doContinue:
            args['--denovo_3dref'] = ''
            args['--sgd_sigma2fudge_halflife'] = self.sgdNoiseVar.get()

    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        if not self.doContinue:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * 2

    def _fillDataFromIter(self, imgSet, iteration):
        outImgsFn = self._getFileName('data', iter=iteration)
        imgSet.setAlignmentProj()
        self.reader = convert.createReader(alignType=ALIGN_PROJ)
        mdIter = convert.Table.iterRows('particles@' + outImgsFn, key='rlnImageId')
        imgSet.copyItems(self._getInputParticles(), doClone=False,
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=mdIter)

    def _createItemMatrix(self, item, row):
        self.reader.setParticleTransform(item, row)

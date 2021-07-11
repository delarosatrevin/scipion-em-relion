# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se)
# *
# * SciLifeLab, Stockholm University
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
import re
from glob import glob
from collections import OrderedDict
from emtable import Table

import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import (BooleanParam, PointerParam, FloatParam,
                                        IntParam, EnumParam, StringParam,
                                        LabelParam, PathParam)
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from pwem.constants import ALIGN_PROJ, ALIGN_NONE
from pwem.emlib.image import ImageHandler
from pwem.objects import SetOfClasses3D, SetOfParticles, SetOfVolumes, Volume
from pwem.protocols import EMProtocol

from relion import Plugin
import relion.convert
from ..constants import ANGULAR_SAMPLING_LIST, MASK_FILL_ZERO


class ProtRelionBase(EMProtocol):
    """ This class contains the common functions for all Relion protocols.
    In subclasses there should be little changes about how to create the command
    line and the files produced.

    Most of the Relion protocols, have two modes: NORMAL or CONTINUE. That's why
    some of the function have a template pattern approach to define the behaviour
    depending on the case.
    """
    _label = None
    IS_CLASSIFY = True
    IS_2D = False
    IS_3D_INIT = False
    OUTPUT_TYPE = SetOfClasses3D
    FILE_KEYS = ['data', 'optimiser', 'sampling']
    CLASS_LABEL = 'rlnClassNumber'
    CHANGE_LABELS = ['rlnChangesOptimalOrientations',
                     'rlnChangesOptimalOffsets',
                     'rlnOverallAccuracyRotations',
                     'rlnOverallAccuracyTranslationsAngst',
                     'rlnChangesOptimalClasses']
    PREFIXES = ['']

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

        self.ClassFnTemplate = '%(rootDir)s/relion_it%(iter)03d_class%(ref)03d.mrc:mrc'
        if not self.doContinue:
            self.continueRun.set(None)
        else:
            if not self.IS_2D:
                self.referenceVolume.set(None)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('relion_it%(iter)03d_')
        myDict = {
            'input_star': self._getPath('input_particles.star'),
            'input_mrcs': self._getPath('input_particles.mrcs'),
            'data_scipion': self.extraIter + 'data_scipion.sqlite',
            'projections': self.extraIter + '%(half)sclass%(ref3d)03d_projections.sqlite',
            'classes_scipion': self.extraIter + 'classes_scipion.sqlite',
            'data': self.extraIter + 'data.star',
            'model': self.extraIter + 'model.star',
            'optimiser': self.extraIter + 'optimiser.star',
            'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
            'all_avgPmax': self._getPath('iterations_avgPmax.star'),
            'all_changes': self._getPath('iterations_changes.star'),
            'selected_volumes': self._getPath('selected_volumes_xmipp.xmd'),
            'dataFinal': self._getExtraPath("relion_data.star"),
            'modelFinal': self._getExtraPath("relion_model.star"),
            'finalvolume': self._getExtraPath("relion_class%(ref3d)03d.mrc:mrc"),
            'final_half1_volume': self._getExtraPath("relion_half1_class%(ref3d)03d_unfil.mrc:mrc"),
            'final_half2_volume': self._getExtraPath("relion_half2_class%(ref3d)03d_unfil.mrc:mrc"),
            'finalSGDvolume': self._getExtraPath("relion_it%(iter)03d_class%(ref3d)03d.mrc:mrc"),
            'preprocess_particles': self._getPath("preprocess_particles.mrcs"),
            'preprocess_particles_star': self._getPath("preprocess_particles.star"),
            'preprocess_particles_prefix': "preprocess_particles"
        }
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FILE_KEYS:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self.PREFIXES:
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc:mrc'

        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000', '???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3})_')

    # -------------------------- DEFINE param functions -----------------------
    def _defineConstants(self):
        self.IS_3D = not self.IS_2D

    def _defineParams(self, form):
        self._defineConstants()

        form.addSection(label='Input')
        # Some hidden variables to be used for conditions
        form.addHidden('isClassify', BooleanParam, default=self.IS_CLASSIFY)
        form.addHidden('is2D', BooleanParam, default=self.IS_2D)

        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                           'run of type *%s* class and most of the input parameters'
                           'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not doContinue',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('copyAlignment', BooleanParam, default=False,
                      label='Consider previous alignment?',
                      condition='not doContinue',
                      help='If set to Yes, then alignment information from'
                           ' input particles will be considered.')
        form.addParam('alignmentAsPriors', BooleanParam, default=False,
                      condition='not doContinue and copyAlignment',
                      expertLevel=LEVEL_ADVANCED,
                      label='Consider alignment as priors?',
                      help='If set to Yes, then alignment information from '
                           'input particles will be considered as PRIORS. This '
                           'option can be used to do restricted local '
                           'search within a range centered around those priors.')
        form.addHidden('fillRandomSubset', BooleanParam, default=True)
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
        form.addParam('maskZero', EnumParam, default=0,
                      choices=['Yes, fill with zeros',
                               'No, fill with random noise'],
                      label='Mask particles with zeros?',
                      condition='not doContinue',
                      help='If set to <Yes>, then in the individual particles, '
                           'the area outside a circle with the radius '
                           'of the particle will be set to zeros prior to '
                           'taking the Fourier transform. '
                           'This will remove noise and therefore increase '
                           'sensitivity in the alignment and classification. '
                           'However, it will also introduce correlations '
                           'between the Fourier components that are not '
                           'modelled. When set to <No>, then the solvent area '
                           'is filled with random noise, which prevents '
                           'introducing correlations.High-resolution '
                           'refinements (e.g. in 3D auto-refine) tend to work '
                           'better when filling the solvent area with random '
                           'noise, some classifications go better when using '
                           'zeros.')
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
        form.addParam('referenceAverages', PointerParam,
                      pointerClass='SetOfAverages', allowsNull=True,
                      condition='not doContinue and isClassify and is2D',
                      expertLevel=LEVEL_ADVANCED,
                      label='Reference averages',
                      help='This option is not recommended and should be used '
                           'with care. The provided averages will be used as '
                           'initial 2D references. If this option is used, '
                           'the number of classes will be ignored. ')

        referenceClass = 'Volume'
        referenceLabel = 'Input volume'
        if self.IS_CLASSIFY:  # allow SetOfVolumes as references for 3D
            referenceClass += ', SetOfVolumes'
            referenceLabel += '(s)'

        if self.IS_3D:
            form.addSection('Reference 3D map',
                            condition='not doContinue and not is2D')

            form.addParam('referenceVolume', PointerParam,
                          pointerClass=referenceClass,
                          important=True,
                          label=referenceLabel,
                          condition='not doContinue and not is2D',
                          help='Initial reference 3D map, it should have the same '
                               'dimensions and the same pixel size as your input '
                               'particles.')

            form.addParam('referenceMask', PointerParam,
                          pointerClass='VolumeMask',
                          label='Reference mask (optional)', allowsNull=True,
                          help='A volume mask containing a (soft) mask with '
                               'the same dimensions as the reference(s), '
                               'and values between 0 and 1, with 1 being 100% '
                               'protein and 0 being 100% solvent. The '
                               'reconstructed reference map will be multiplied '
                               'by this mask. If no mask is given, a soft '
                               'spherical mask based on the <radius> of the '
                               'mask for the experimental images will be '
                               'applied.\n\n'
                               'In some cases, for example for non-empty '
                               'icosahedral viruses, it is also useful to use '
                               'a second mask. Check _Advaced_ level and '
                               'select another volume mask')
            form.addParam('solventMask', PointerParam, pointerClass='VolumeMask',
                          expertLevel=LEVEL_ADVANCED, allowsNull=True,
                          label='Second reference mask (optional)',
                          help='For all white (value 1) pixels in this second '
                               'mask the corresponding pixels in the '
                               'reconstructed map are set to the average value '
                               'of these pixels. Thereby, for example, the '
                               'higher density inside the virion may be set to '
                               'a constant. Note that this second mask should '
                               'have one-values inside the virion and '
                               'zero-values in the capsid and the solvent '
                               'areas.')
            form.addParam('solventFscMask', BooleanParam, default=False,
                          condition='not isClassify',
                          label='Use solvent-flattened FSCs?',
                          help='If set to Yes, then instead of using '
                               'unmasked maps to calculate the gold-standard '
                               'FSCs during refinement, masked half-maps '
                               'are used and a post-processing-like '
                               'correction of the FSC curves (with '
                               'phase-randomisation) is performed every '
                               'iteration. This only works when a reference '
                               'mask is provided on the I/O tab. This may '
                               'yield higher-resolution maps, especially '
                               'when the mask contains only a relatively '
                               'small volume inside the box.')
            form.addParam('isMapAbsoluteGreyScale', BooleanParam, default=False,
                          condition='not doContinue',
                          label="Is initial 3D map on absolute greyscale?",
                          help='The probabilities are based on squared differences,'
                               ' so that the absolute grey scale is important. \n'
                               'Probabilities are calculated based on a Gaussian '
                               'noise model, which contains a squared difference '
                               'term between the reference and the experimental '
                               'image. This has a consequence that the reference '
                               'needs to be on the same absolute intensity '
                               'grey-scale as the experimental images. RELION and '
                               'XMIPP reconstruct maps at their absolute '
                               'intensity grey-scale. Other packages may perform '
                               'internal normalisations of the reference density, '
                               'which will result in incorrect grey-scales. '
                               'Therefore: if the map was reconstructed in RELION '
                               'or in XMIPP, set this option to Yes, otherwise '
                               'set it to No. If set to No, RELION will use a ('
                               'grey-scale invariant) cross-correlation criterion '
                               'in the first iteration, and prior to the second '
                               'iteration the map will be filtered again using '
                               'the initial low-pass filter. This procedure is '
                               'relatively quick and typically does not '
                               'negatively affect the outcome of the subsequent '
                               'MAP refinement. Therefore, if in doubt it is '
                               'recommended to set this option to No.')

            self.addSymmetry(form)

            form.addParam('initialLowPassFilterA', FloatParam, default=60,
                          condition='not is2D and not doContinue',
                          label='Initial low-pass filter (A)',
                          help='It is recommended to strongly low-pass filter your '
                               'initial reference map. If it has not yet been '
                               'low-pass filtered, it may be done internally using '
                               'this option. If set to 0, no low-pass filter will '
                               'be applied to the initial reference(s).')
        else:
            form.addParam('referenceMask2D', PointerParam, pointerClass='Mask',
                          label='Reference mask (optional)', allowsNull=True,
                          expertLevel=LEVEL_ADVANCED,
                          help='User-provided mask for the references ('
                               'default is to use spherical mask with '
                               'particle_diameter)')

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
        if not Plugin.IS_GT31():
            form.addParam('hasReferenceCTFCorrected', BooleanParam, default=False,
                          condition='not is2D and not doContinue',
                          label='Has reference been CTF-corrected?',
                          help='Set this option to Yes if the reference map '
                               'represents CTF-unaffected density, e.g. it was '
                               'created using Wiener filtering inside RELION or '
                               'from a PDB. If set to No, then in the first '
                               'iteration, the Fourier transforms of the reference '
                               'projections are not multiplied by the CTFs.')
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

        if self.IS_CLASSIFY:
            # In our organization of the parameters, Optimisation tab only
            # make sense when in a classification case
            form.addSection(label='Optimisation')

            form.addParam('numberOfClasses', IntParam, default=3,
                          condition='not doContinue and isClassify',
                          label='Number of classes:',
                          help='The number of classes (K) for a multi-reference '
                               'refinement. These classes will be made in an '
                               'unsupervised manner from a single reference by '
                               'division of the data into random subsets during '
                               'the first iteration.')
            # Default T is 2 for 2D but 4 for 3D in Relion GUI
            form.addParam('regularisationParamT', FloatParam,
                          default=2 if self.IS_2D else 4,
                          label='Regularisation parameter T',
                          help='Bayes law strictly determines the relative '
                               'weight between the contribution of the '
                               'experimental data and the prior. '
                               'However, in practice one may need to adjust '
                               'this weight to put slightly more weight on the '
                               'experimental data to allow optimal results. '
                               'Values greater than 1 for this regularisation '
                               'parameter (T in the JMB2011 paper) put more '
                               'weight on the experimental data. Values around '
                               '2-4 have been observed to be useful for 3D '
                               'refinements, values of 1-2 for 2D refinements. '
                               'Too small values yield too-low resolution '
                               'structures; too high values result in '
                               'over-estimated resolutions and overfitting.')
            form.addParam('numberOfIterations', IntParam, default=25,
                          label='Number of iterations',
                          help='Number of iterations to be performed. Note '
                               'that the current implementation does NOT '
                               'comprise a convergence criterium. Therefore, '
                               'the calculations will need to be stopped '
                               'by the user if further iterations do not yield '
                               'improvements in resolution or classes. '
                               'If continue option is True, you going to do '
                               'this number of new iterations (e.g. if '
                               '*Continue from iteration* is set 3 and this '
                               'param is set 25, the final iteration of the '
                               'protocol will be the 28th.')

            if (Plugin.IS_GT31() and self.IS_3D) or not Plugin.IS_GT31():
                form.addParam('useFastSubsets', BooleanParam, default=False,
                              condition='not doContinue',
                              label='Use fast subsets (for large data sets)?',
                              help='If set to Yes, the first 5 iterations will '
                                   'be done with random subsets of only K*100 '
                                   'particles (K being the number of classes); '
                                   'the next 5 with K*300 particles, the next '
                                   '5 with 30% of the data set; and the final '
                                   'ones with all data. This was inspired by '
                                   'a cisTEM implementation by Niko Grigorieff'
                                   ' et al.')

            if Plugin.IS_GT31() and self.IS_2D:  # relion4 2D cls case
                form.addParam('useGradientAlg', BooleanParam, default=False,
                              condition='not doContinue',
                              label='Use gradient-driven algorithm?',
                              help='If set to Yes, use faster NGrad '
                                   'algorithm instead of the default '
                                   'expectation maximization. If used, *increase '
                                   'number of iterations to about 100!*')
                form.addParam('centerAvg', BooleanParam, default=True,
                              label='Center class averages?',
                              help='If set to Yes, every iteration the class '
                                   'average images will be centered on their '
                                   'center-of-mass. This will work only for '
                                   'positive signals, so the particles should '
                                   'be white.')

            form.addParam('limitResolEStep', FloatParam, default=-1,
                          label='Limit resolution E-step to (A)',
                          condition="not doContinue",
                          help='If set to a positive number, then the '
                               'expectation step (i.e. the alignment) will be '
                               'done only including the Fourier components up '
                               'to this resolution (in Angstroms). This is '
                               'useful to prevent overfitting, as the '
                               'classification runs in RELION are not to be '
                               'guaranteed to be 100% overfitting-free (unlike '
                               'the _3D auto-refine_ with its gold-standard '
                               'FSC). In particular for very difficult data '
                               'sets, e.g. of very small or featureless '
                               'particles, this has been shown to give much '
                               'better class averages. In such cases, values '
                               'in the range of 7-12 Angstroms have proven '
                               'useful.')

            # Change the Sampling section name depending if classify or refine 3D
            form.addSection('Sampling')
        else:
            form.addSection('Auto-Sampling')

        form.addParam('doImageAlignment', BooleanParam, default=True,
                      label='Perform image alignment?', condition="isClassify",
                      help='If set to No, then rather than performing both alignment '
                           'and classification, only classification will be performed. '
                           'This allows the use of very focused masks. This requires '
                           'that the optimal orientations of all particles are already '
                           'calculated.')
        if self.IS_3D:
            form.addParam('angularSamplingDeg', EnumParam, default=2,
                          choices=ANGULAR_SAMPLING_LIST,
                          label='Initial angular sampling (deg)',
                          condition='doImageAlignment and (isClassify or not doContinue)',
                          help='There are only a few discrete angular samplings'
                               ' possible because we use the HealPix library to'
                               ' generate the sampling of the first two Euler '
                               'angles on the sphere. The samplings are '
                               'approximate numbers and vary slightly over '
                               'the sphere.')
        else:
            form.addParam('inplaneAngularSamplingDeg', FloatParam, default=6,
                          label='In-plane angular sampling (deg)',
                          condition="doImageAlignment",
                          help='The sampling rate for the in-plane rotation '
                               'angle (psi) in degrees.\n'
                               'Using fine values will slow down the program. '
                               'Recommended value for\n'
                               'most 2D refinements: 5 degrees. \n\n'
                               'If auto-sampling is used, this will be the '
                               'value for the first \niteration(s) only, and '
                               'the sampling rate will be increased \n'
                               'automatically after that.')
        form.addParam('offsetSearchRangePix', FloatParam, default=5,
                      condition='doImageAlignment and (isClassify or not doContinue)',
                      label='Initial offset range (pix)',
                      help='Probabilities will be calculated only for '
                           'translations in a circle with this radius (in '
                           'pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal '
                           'translation for each image in the previous '
                           'iteration.')
        form.addParam('offsetSearchStepPix', FloatParam, default=1.0,
                      condition='doImageAlignment and (isClassify or not doContinue)',
                      label='Initial offset step (pix)',
                      help='Translations will be sampled with this step-size '
                           '(in pixels). Translational sampling is also done '
                           'using the adaptive approach. Therefore, if '
                           'adaptive=1, the translations will first be '
                           'evaluated on a 2x coarser grid.')
        if self.IS_3D:
            if self.IS_CLASSIFY:
                form.addParam('localAngularSearch', BooleanParam, default=False,
                              condition='not is2D and doImageAlignment',
                              label='Perform local angular search?',
                              help='If set to Yes, then rather than performing '
                                   'exhaustive angular searches, local '
                                   'searches within the range given below will '
                                   'be performed. A prior Gaussian distribution'
                                   ' centered at the optimal orientation in the'
                                   ' previous iteration and with a stddev of '
                                   '1/3 of the range given below will be '
                                   'enforced.')
                form.addParam('localAngularSearchRange', FloatParam,
                              default=5.0,
                              condition='localAngularSearch and doImageAlignment',
                              label='Local angular search range',
                              help='Local angular searches will be performed '
                                   'within +/- the given amount (in degrees) '
                                   'from the optimal orientation in the '
                                   'previous iteration. A Gaussian prior (also '
                                   'see previous option) will be applied, so '
                                   'that orientations closer to the optimal '
                                   'orientation in the previous iteration will '
                                   'get higher weights than those further away.')
                form.addParam('relaxSymm', StringParam, default='',
                              condition='localAngularSearch and doImageAlignment',
                              label='Relax symmetry',
                              help="With this option, poses related to the standard "
                                   "local angular search range by the given point "
                                   "group will also be explored. For example, if "
                                   "you have a pseudo-symmetric dimer A-A', "
                                   "refinement or classification in C1 with symmetry "
                                   "relaxation by C2 might be able to improve "
                                   "distinction between A and A'. Note that the "
                                   "reference must be more-or-less aligned to the "
                                   "convention of (pseudo-)symmetry operators. "
                                   "For details, see Ilca et al 2019 and "
                                   "Abrishami et al 2020.")
            else:
                form.addParam('localSearchAutoSamplingDeg', EnumParam,
                              condition='not doContinue',
                              default=4, choices=ANGULAR_SAMPLING_LIST,
                              label='Local search from auto-sampling (deg)',
                              help='In the automated procedure to increase the '
                                   'angular samplings, local angular searches '
                                   'of -6/+6 times the sampling rate will be '
                                   'used from this angular sampling rate '
                                   'onwards.')
                form.addParam('relaxSymm', StringParam, default='',
                              condition='doImageAlignment',
                              label='Relax symmetry',
                              help="With this option, poses related to the standard "
                                   "local angular search range by the given point "
                                   "group will also be explored. For example, if "
                                   "you have a pseudo-symmetric dimer A-A', "
                                   "refinement or classification in C1 with symmetry "
                                   "relaxation by C2 might be able to improve "
                                   "distinction between A and A'. Note that the "
                                   "reference must be more-or-less aligned to the "
                                   "convention of (pseudo-)symmetry operators. "
                                   "For details, see Ilca et al 2019 and "
                                   "Abrishami et al 2020.")
                form.addParam('useFinerSamplingFaster', BooleanParam,
                              default=False,
                              label='Use finer angular sampling faster?',
                              help='If set to Yes, then let auto-refinement '
                                   'proceed faster with finer angular '
                                   'samplings. Two additional command-line '
                                   'options will be passed to the refine '
                                   'program:\n\n'
                                   '\t--auto_ignore_angles lets angular '
                                   'sampling go down despite changes '
                                   'still happening in the angles\n'
                                   '\t--auto_resol_angles lets angular '
                                   'sampling go down if the current '
                                   'resolution already requires that '
                                   'sampling at the edge of the particle.\n\n'
                                   'This option will make the computation '
                                   'faster, but has not been tested for '
                                   'many cases for potential loss in '
                                   'reconstruction quality upon convergence.')



        if self.IS_CLASSIFY:
            form.addParam('allowCoarserSampling', BooleanParam,
                          condition='doImageAlignment',
                          default=False,
                          label='Allow coarser sampling?',
                          help='If set to Yes, the program will use '
                               'coarser angular and translational '
                               'samplings if the estimated accuracies '
                               'of the assignments is still low in the '
                               'earlier iterations. This may speed up '
                               'the calculations.')

        form.addSection('Compute')
        self._defineComputeParams(form)

        joinHalves = ("--low_resol_join_halves 40 (only not continue mode)"
                      if not self.IS_CLASSIFY else "")

        form.addParam('oversampling', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Over-sampling",
                      help="Adaptive oversampling order to speed-up "
                           "calculations (0=no oversampling, 1=2x, 2=4x, etc)")

        form.addParam('extraParams', StringParam,
                      default='',
                      label='Additional arguments',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program, e.g: \n"
                           "--dont_combine_weights_via_disc\n"
                           "--verb 1\n"
                           "--pad 2\n" + joinHalves)

        form.addParallelSection(threads=1, mpi=3)

    def addSymmetry(self, container):
        container.addParam('symmetryGroup', StringParam, default='c1',
                           condition='not doContinue',
                           label="Symmetry",
                           help='See [[https://relion.readthedocs.io/'
                                'en/latest/Reference/Conventions.html#symmetry]'
                                '[Relion Symmetry]] '
                                'page for a description of the symmetry format '
                                'accepted by Relion')

    def _defineComputeParams(self, form):
        form.addParam('useParallelDisk', BooleanParam, default=True,
                      label='Use parallel disc I/O?',
                      help='If set to Yes, all MPI slaves will read '
                           'their own images from disc. Otherwise, only '
                           'the master will read images and send them '
                           'through the network to the slaves. Parallel '
                           'file systems like gluster of fhgfs are good '
                           'at parallel disc I/O. NFS may break with many '
                           'slaves reading in parallel.')
        form.addParam('pooledParticles', IntParam, default=3,
                      label='Number of pooled particles:',
                      help='Particles are processed in individual batches '
                           'by MPI slaves. During each batch, a stack of '
                           'particle images is only opened and closed '
                           'once to improve disk access times. All '
                           'particle images of a single batch are read '
                           'into memory together. The size of these '
                           'batches is at least one particle per thread '
                           'used. The nr_pooled_particles parameter '
                           'controls how many particles are read together '
                           'for each thread. If it is set to 3 and one '
                           'uses 8 threads, batches of 3x8=24 particles '
                           'will be read together. This may improve '
                           'performance on systems where disk access, and '
                           'particularly metadata handling of disk '
                           'access, is a problem. It has a modest cost of '
                           'increased RAM usage.')
        if self.IS_3D:
            if not (Plugin.IS_GT31() and self.IS_3D_INIT):
                form.addParam('skipPadding', BooleanParam, default=False,
                              label='Skip padding',
                              help='If set to Yes, the calculations will not use '
                                   'padding in Fourier space for better '
                                   'interpolation in the references. Otherwise, '
                                   'references are padded 2x before Fourier '
                                   'transforms are calculated. Skipping padding '
                                   '(i.e. use --pad 1) gives nearly as good results '
                                   'as using --pad 2, but some artifacts may appear '
                                   'in the corners from signal that is folded back.')
            form.addParam('skipGridding', BooleanParam, default=True,
                          label='Skip gridding',
                          help='If set to Yes, the calculations will '
                               'skip gridding in the M step to save time, '
                               'typically with just as good results.')

        form.addParam('allParticlesRam', BooleanParam, default=False,
                      label='Pre-read all particles into RAM?',
                      help='If set to Yes, all particle images will be '
                           'read into computer memory, which will greatly '
                           'speed up calculations on systems with slow '
                           'disk access. However, one should of course be '
                           'careful with the amount of RAM available. '
                           'Because particles are read in '
                           'float-precision, it will take \n'
                           '( N * (box_size)^2 * 4 / (1024 * 1024 '
                           '* 1024) ) Giga-bytes to read N particles into '
                           'RAM. For 100 thousand 200x200 images, that '
                           'becomes 15Gb, or 60 Gb for the same number of '
                           '400x400 particles. Remember that running a '
                           'single MPI slave on each node that runs as '
                           'many threads as available cores will have '
                           'access to all available RAM.\n\n'
                           'If parallel disc I/O is set to No, then only '
                           'the master reads all particles into RAM and '
                           'sends those particles through the network to '
                           'the MPI slaves during the refinement '
                           'iterations.')
        form.addParam('scratchDir', PathParam,
                      condition='not allParticlesRam',
                      label='Copy particles to scratch directory: ',
                      help='If a directory is provided here, then the job '
                           'will create a sub-directory in it called '
                           'relion_volatile. If that relion_volatile '
                           'directory already exists, it will be wiped. '
                           'Then, the program will copy all input '
                           'particles into a large stack inside the '
                           'relion_volatile subdirectory. Provided this '
                           'directory is on a fast local drive (e.g. an '
                           'SSD drive), processing in all the iterations '
                           'will be faster. If the job finishes '
                           'correctly, the relion_volatile directory will '
                           'be wiped. If the job crashes, you may want to '
                           'remove it yourself.')
        form.addParam('combineItersDisc', BooleanParam, default=False,
                      label='Combine iterations through disc?',
                      help='If set to Yes, at the end of every iteration '
                           'all MPI slaves will write out a large file '
                           'with their accumulated results. The MPI '
                           'master will read in all these files, combine '
                           'them all, and write out a new file with the '
                           'combined results. All MPI slaves will then '
                           'read in the combined results. This reduces '
                           'heavy load on the network, but increases load '
                           'on the disc I/O. This will affect the time it '
                           'takes between the progress-bar in the '
                           'expectation step reaching its end (the mouse '
                           'gets to the cheese) and the start of the '
                           'ensuing maximisation step. It will depend on '
                           'your system setup which is most efficient.')
        form.addParam('doGpu', BooleanParam, default=True,
                      label='Use GPU acceleration?',
                      help='If set to Yes, the job will try to use GPU '
                           'acceleration.')
        form.addParam('gpusToUse', StringParam, default='',
                      label='Which GPUs to use:', condition='doGpu',
                      help='This argument is not necessary. If left empty, '
                           'the job itself will try to allocate available '
                           'GPU resources. You can override the default '
                           'allocation by providing a list of which GPUs '
                           '(0,1,2,3, etc) to use. MPI-processes are '
                           'separated by ":", threads by ",". '
                           'For example: "0,0:1,1:0,0:1,1"')

    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('convertInputStep',
                                 self._getInputParticles().getObjId(),
                                 bool(self.copyAlignment))
        self._insertRelionStep()
        self._insertFunctionStep('createOutputStep')

    def _insertRelionStep(self):
        """ Prepare the command line arguments before calling Relion. """
        # Join in a single line all key, value pairs of the args dict
        args = OrderedDict()

        if self.doContinue:
            self._setContinueArgs(args)
        else:
            self._setNormalArgs(args)

        self._setComputeArgs(args)

        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])

        if self.extraParams.hasValue():
            params += ' ' + self.extraParams.get()

        self._insertFunctionStep('runRelionStep', params)

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, particlesId, copyAlignment):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId: use this parameters just to force redo of convert if
                the input particles are changed.
        """
        imgSet = self._getInputParticles()
        if not self.doContinue:
            imgStar = self._getFileName('input_star')

            self.info("Converting set from '%s' into '%s'" %
                      (imgSet.getFileName(), imgStar))

            # Pass stack file as None to avoid write the images files
            # If copyAlignment is set to False pass alignType to ALIGN_NONE
            alignType = imgSet.getAlignment() if copyAlignment else ALIGN_NONE
            hasAlign = alignType != ALIGN_NONE
            alignToPrior = hasAlign and getattr(self, 'alignmentAsPriors', False)

            if self.doCtfManualGroups:
                self._defocusGroups = self.createDefocusGroups()
                print(self._defocusGroups)

            relion.convert.writeSetOfParticles(
                imgSet, imgStar,
                outputDir=self._getExtraPath(),
                alignType=alignType,
                postprocessImageRow=self._postprocessParticleRow)

            if alignToPrior:
                mdOptics = Table(fileName=imgStar, tableName='optics')
                mdParts = Table(fileName=imgStar, tableName='particles')
                self._copyAlignAsPriors(mdParts, alignType)

                with open(imgStar, "w") as f:
                    mdParts.writeStar(f, tableName='particles')
                    mdOptics.writeStar(f, tableName='optics')

            if self._getRefArg():
                self._convertRef()
        else:
            self.info("In continue mode is not necessary convert the input "
                      "particles")

    def runRelionStep(self, params):
        """ Execute the relion steps with the give params. """
        params += ' --j %d' % self.numberOfThreads
        self.runJob(self._getProgram(), params)

    def _getEnviron(self):
        env = Plugin.getEnviron()

        if self.useGpu():
            prepend = env.get('RELION_PREPEND', '')
        else:
            prepend = env.get('RELION_PREPEND_CPU', '')

        env.setPrepend(prepend)

        return env

    def createOutputStep(self):
        pass  # should be implemented in subclasses

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        if self.doContinue:
            continueProtocol = self.continueRun.get()
            if (continueProtocol is not None and
                    continueProtocol.getObjId() == self.getObjId()):
                errors.append('In Scipion you must create a new Relion run')
                errors.append('and select the continue option rather than')
                errors.append('select continue from the same run.')
                errors.append('')  # add a new line
            errors += self._validateContinue()
        else:
            if self._getInputParticles().isOddX():
                errors.append("Relion only works with even values for the "
                              "image dimensions!")

            # if doing scaling, the input is not on abs greyscale
            if self.getAttributeValue('referenceVolume'):
                volX = self._getReferenceVolumes()[0].getXDim()
                ptclX = self._getInputParticles().getXDim()
                if (ptclX != volX) and self.isMapAbsoluteGreyScale:
                    errors.append("Input particles and references have "
                                  "different dimensions, so the reference is "
                                  "not on the absolute greyscale. Select *No* for "
                                  "that option to continue.")

            errors += self._validateNormal()

        if self.IS_CLASSIFY and not self.doImageAlignment:
            if self.doGpu:
                errors.append('When only doing classification (no alignment) '
                              'GPU acceleration can not be used.')
            if not self.copyAlignment:
                errors.append('If you want to do only classification without '
                              'alignment, then you should use the option: \n'
                              '*Consider previous alignment?* = Yes')
            if self.alignmentAsPriors:
                errors.append('When only doing classification (no alignment) '
                              " option *Consider alignment as priors* cannot"
                              " be enabled.")

        return errors

    def _validateNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _validateContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []

    def _citations(self):
        cites = []
        return cites

    def _summary(self):
        self._initialize()

        lastIter = self._lastIter()

        if lastIter is not None:
            iterMsg = 'Iteration %d' % lastIter
            if self.hasAttribute('numberOfIterations'):
                iterMsg += '/%d' % self._getnumberOfIters()
        else:
            iterMsg = 'No iteration finished yet.'
        summary = [iterMsg]

        inputParts = self._getInputParticles()
        if inputParts is not None and inputParts.isPhaseFlipped():
            summary.append("Your input images are ctf-phase flipped")

        if self.doContinue:
            summary += self._summaryContinue()
        summary += self._summaryNormal()
        return summary

    def _summaryNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _summaryContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []

    def _methods(self):
        """ Should be overwritten in each protocol.
        """
        return []

    # -------------------------- UTILS functions ------------------------------
    def _setNormalArgs(self, args):
        inputParts = self._getInputParticles()
        ps = inputParts.getSamplingRate()

        maskDiameter = self.maskDiameterA.get()
        if maskDiameter <= 0:
            maskDiameter = ps * inputParts.getXDim()

        args['--i'] = self._getFileName('input_star')
        args['--particle_diameter'] = maskDiameter

        self._setCTFArgs(args)

        if self.maskZero == MASK_FILL_ZERO:
            args['--zero_mask'] = ''

        if self.IS_CLASSIFY:
            args['--K'] = self.numberOfClasses.get()
            if self.limitResolEStep > 0:
                args['--strict_highres_exp'] = self.limitResolEStep.get()
            if self.IS_2D and Plugin.IS_GT31():
                if self.useGradientAlg:
                    args['--grad'] = ''
                    args['--init_blobs'] = ''
                    args['--class_inactivity_threshold'] = 0.1
                if self.centerAvg:
                    args['--center_classes'] = ''

        if self.IS_3D and not self.IS_3D_INIT:
            if not self.isMapAbsoluteGreyScale:
                args['--firstiter_cc'] = ''
            args['--ini_high'] = self.initialLowPassFilterA.get()
            args['--sym'] = self.symmetryGroup.get()
            # We use the same pixel size as input particles, since
            # we convert anyway the input volume to match same size
            args['--ref_angpix'] = ps

        refArg = self._getRefArg()
        if refArg:
            args['--ref'] = refArg

        self._setBasicArgs(args)

    def _setContinueArgs(self, args):
        continueRun = self.continueRun.get()
        continueRun._initialize()

        self._setBasicArgs(args)

        continueIter = self._getContinueIter()
        args['--continue'] = continueRun._getFileName('optimiser',
                                                      iter=continueIter)

    def _getScratchDir(self):
        """ Returns the scratch dir value without spaces.
         If none, the empty string will be returned.
        """
        scratchDir = self.scratchDir.get() or ''
        return scratchDir.strip()

    def _setComputeArgs(self, args):
        if not self.combineItersDisc:
            args['--dont_combine_weights_via_disc'] = ''

        if not self.useParallelDisk:
            args['--no_parallel_disc_io'] = ''

        if self.allParticlesRam:
            args['--preread_images'] = ''
        else:
            if self._getScratchDir():
                args['--scratch_dir'] = self._getScratchDir()

        args['--pool'] = self.pooledParticles.get()

        if self.doGpu:
            gpuStr = self.gpusToUse.get().strip()
            if gpuStr:
                if not gpuStr.startswith('"'):
                    gpuStr = '"' + gpuStr
                if not gpuStr.endswith('"'):
                    gpuStr += '"'

            args['--gpu'] = gpuStr

    def _getSamplingFactor(self):
        return 1 if self.oversampling == 0 else 2 * self.oversampling.get()

    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args['--norm'] = ''
        args['--scale'] = ''
        args['--o'] = self._getExtraPath('relion')
        args['--oversampling'] = self.oversampling.get()
        args['--flatten_solvent'] = ''

        if self.IS_CLASSIFY:
            args['--tau2_fudge'] = self.regularisationParamT.get()
            args['--iter'] = self._getnumberOfIters()

            if not self.doContinue:
                self._setSubsetArgs(args)

        # Padding can be set in a normal run or in a continue
        if self.IS_3D:
            args['--pad'] = 1 if self.skipPadding else 2
            if self.skipGridding:
                args['--skip_gridding'] = ''

        self._setSamplingArgs(args)
        self._setMaskArgs(args)

    def _setSamplingArgs(self, args):
        """ Should be implemented in subclasses. """
        pass

    def _setCTFArgs(self, args):
        # CTF stuff
        if self.doCTF:
            args['--ctf'] = ''

            # this only can be true if is 3D
            if not Plugin.IS_GT31() and self.hasReferenceCTFCorrected:
                args['--ctf_corrected_ref'] = ''

            if self._getInputParticles().isPhaseFlipped():
                args['--ctf_phase_flipped'] = ''

            if self.ignoreCTFUntilFirstPeak:
                args['--ctf_intact_first_peak'] = ''

    def _setMaskArgs(self, args):
        if self.IS_3D:
            tmp = self._getTmpPath()
            newDim = self._getInputParticles().getXDim()
            newPix = self._getInputParticles().getSamplingRate()
            if self.referenceMask.hasValue():
                mask = relion.convert.convertMask(self.referenceMask.get(),
                                                  tmp, newPix, newDim)
                args['--solvent_mask'] = mask

            if self.solventMask.hasValue():
                solventMask = relion.convert.convertMask(self.solventMask.get(),
                                                         tmp, newPix, newDim)
                args['--solvent_mask2'] = solventMask

            if self.referenceMask.hasValue() and self.solventFscMask:
                args['--solvent_correct_fsc'] = ''
        else:
            if self.referenceMask2D.hasValue():
                tmp = self._getTmpPath()
                newDim = self._getInputParticles().getXDim()
                newPix = self._getInputParticles().getSamplingRate()
                mask = relion.convert.convertMask(self.referenceMask2D.get(),
                                                  tmp, newPix, newDim)
                args['--solvent_mask'] = mask

    def _setSubsetArgs(self, args):
        if self._useFastSubsets():
            args['--fast_subsets'] = ''

    def _getProgram(self, program='relion_refine'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'

        return program

    def _runProgram(self, program, args, **kwargs):
        """ Helper function to get the program name if mpi are used and
        call runJob function.
        """
        return self.runJob(self._getProgram(program), args, **kwargs)

    def _getInputParticles(self):
        if self.doContinue:
            self.inputParticles.set(self.continueRun.get().inputParticles.get())
        return self.inputParticles.get()

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 3 digits iteration number
        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)

        if clean:
            pwutils.cleanPath(data_classes)

        if not os.path.exists(data_classes):
            clsSet = self.OUTPUT_TYPE(filename=data_classes)
            clsSet.setImages(self.inputParticles.get())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _fillClassesFromIter(self, clsSet, iteration):
        """ Should be implemented in subclasses. """
        pass

    def _getIterData(self, it, **kwargs):
        """ Sort the it??.data.star file by the maximum likelihood. """
        data_sqlite = self._getFileName('data_scipion', iter=it)

        if not os.path.exists(data_sqlite):
            iterImgSet = SetOfParticles(filename=data_sqlite)
            iterImgSet.copyInfo(self._getInputParticles())
            self._fillDataFromIter(iterImgSet, it)
            iterImgSet.write()
            iterImgSet.close()

        return data_sqlite

    def _fillDataFromIter(self, imgSet, iteration):
        """ Should be implemented in subclasses. """
        pass

    def _getContinueIter(self):
        continueRun = self.continueRun.get()

        if continueRun is not None:
            continueRun._initialize()

        if self.doContinue:
            if self.continueIter.get() == 'last':
                continueIter = continueRun._lastIter()
            else:
                continueIter = int(self.continueIter.get())
        else:
            continueIter = 0

        return continueIter

    def _getnumberOfIters(self):
        return self._getContinueIter() + self.numberOfIterations.get()

    def _getReferenceVolumes(self):
        """ Return a list with all input references.
        (Could be one or more volumes).
        """
        inputObj = self.referenceVolume.get()

        if isinstance(inputObj, Volume):
            return [inputObj]
        elif isinstance(inputObj, SetOfVolumes):
            return [vol.clone() for vol in inputObj]
        else:
            raise Exception("Invalid input reference of class: %s"
                            % inputObj.getClassName())

    def _getRefArg(self):
        """ Return the filename that will be used for the --ref argument.
        The value will depend if in 2D and 3D or if input references will
        be used.
        It will return None if no --ref should be used. """
        if self.IS_3D:
            if not self.IS_3D_INIT:
                refVols = self._getReferenceVolumes()
                if len(refVols) == 1:
                    return self._convertVolFn(refVols[0])
                else:  # input SetOfVolumes as references
                    return self._getRefStar()
        else:  # 2D
            if self.referenceAverages.get():
                return self._getRefStar()
        return None  # No --ref should be used at this point

    def _convertVolFn(self, inputVol):
        """ Return a new name if the inputFn is not .mrc """
        index, fn = inputVol.getLocation()
        return self._getTmpPath(pwutils.replaceBaseExt(fn, '%02d.mrc' % index))

    def _convertVol(self, inputVol):
        outputFn = self._convertVolFn(inputVol)

        if outputFn:
            newPix = self._getInputParticles().getSamplingRate()
            newDim = self._getInputParticles().getXDim()
            relion.convert.convertMask(inputVol, outputFn, newPix=newPix,
                                       newDim=newDim, threshold=False)

        return outputFn

    def _getRefStar(self):
        return self._getTmpPath("input_references.star")

    def _convertRef(self):
        ih = ImageHandler()

        if self.IS_3D:
            if not self.IS_3D_INIT:
                refVols = self._getReferenceVolumes()
                if len(refVols) == 1:
                    self._convertVol(refVols[0])
                else:  # input SetOfVolumes as references
                    table = Table(columns=['rlnReferenceImage'])
                    for vol in refVols:
                        newVolFn = self._convertVol(vol)
                        table.addRow(newVolFn)
                    with open(self._getRefStar(), 'w') as f:
                        table.writeStar(f)
        else:  # 2D
            inputAvgs = self.referenceAverages.get()
            if inputAvgs:
                table = Table(columns=['rlnReferenceImage'])
                refStack = self._getTmpPath('input_references.mrcs')
                for i, avg in enumerate(inputAvgs):
                    newAvgLoc = (i + 1, refStack)
                    ih.convert(avg, newAvgLoc)
                    table.addRow("%05d@%s" % newAvgLoc)
                with open(self._getRefStar(), 'w') as f:
                    table.writeStar(f)

    def _postprocessParticleRow(self, part, partRow):
        if self.doCtfManualGroups:
            groupId = self._defocusGroups.getGroup(part.getCTF().getDefocusU()).id
            partRow['rlnGroupName'] = "ctf_group_%03d" % groupId

    def _useFastSubsets(self):
        return self.getAttributeValue('useFastSubsets', False)

    def useGpu(self):
        """
        Return True if the protocol has gpu option and it has been selected.
        """
        return self.getAttributeValue('doGpu', False)

    def _copyAlignAsPriors(self, mdParts, alignType):
        # set priors equal to orig. values
        mdParts.addColumns('rlnOriginXPriorAngst=rlnOriginXAngst')
        mdParts.addColumns('rlnOriginYPriorAngst=rlnOriginYAngst')
        mdParts.addColumns('rlnAnglePsiPrior=rlnAnglePsi')

        if alignType == ALIGN_PROJ:
            mdParts.addColumns('rlnAngleRotPrior=rlnAngleRot')
            mdParts.addColumns('rlnAngleTiltPrior=rlnAngleTilt')

    def createDefocusGroups(self):
        defocusGroups = None

        if self.doCtfManualGroups:
            defocusGroups = relion.convert.DefocusGroups()
            defocusGroups.splitByDiff(self._getInputParticles(),
                                      defocusDiff=self.defocusRange.get(),
                                      minGroupSize=self.numParticles.get())
        return defocusGroups

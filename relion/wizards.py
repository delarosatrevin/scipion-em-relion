# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os

from pwem.constants import UNIT_PIXEL, UNIT_ANGSTROM, FILTER_LOW_PASS_NO_DECAY
from pwem.viewers import EmPlotter
from pwem.wizards.wizard import (ParticleMaskRadiusWizard, FilterVolumesWizard,
                                 EmWizard, ColorScaleWizardBase,
                                 BandPassFilterDialog, dialog)
from pyworkflow.gui.browser import FileBrowserWindow

import relion.convert as convert
from .protocols import *
from .viewers import RelionLocalResViewer

# =============================================================================
# MASKS
# =============================================================================


class RelionBackRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(ProtRelionPreprocessParticles, ['backRadius'])]
    _unit = UNIT_PIXEL
    
    def _getProtocolImages(self, protocol):
        return protocol.inputParticles
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        return {
            'input': self._getProtocolImages(protocol),
            'label': label,
            'value': value
        }

    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form, *args):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, units=self._unit)

    def setVar(self, form, label, value):
        form.setVar(label, int(value))


class RelionPartMaskDiameterWizard(RelionBackRadiusWizard):
    _targets = [(ProtRelionClassify2D, ['maskDiameterA']),
                (ProtRelionRefine3D, ['maskDiameterA']),
                (ProtRelionClassify3D, ['maskDiameterA']),
                (ProtRelionClassify2D, ['maskDiameterA']),
                (ProtRelionInitialModel, ['maskDiameterA'])]
    _unit = UNIT_ANGSTROM

    def _getParameters(self, protocol):
        protParams = RelionBackRadiusWizard._getParameters(self, protocol)
        # adjust to from diameter to radius
        protParams['value'] /= 2

        return protParams

    def setVar(self, form, label, value):
        # adjust again from radius to diameter
        form.setVar(label, value * 2)


# =============================================================================
# FILTER
# =============================================================================

class RelionVolFilterWizard(FilterVolumesWizard):
    _targets = [(ProtRelionClassify3D, ['initialLowPassFilterA']),
                (ProtRelionRefine3D, ['initialLowPassFilterA']),
                (ProtRelionCreateMask3D, ['initialLowPassFilterA'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        protParams = {}

        if protocol.__class__.__name__ == 'ProtRelionCreateMask3D':
            protParams['input'] = protocol.inputVolume
        else:
            protParams['input'] = protocol.referenceVolume

        protParams['label'] = label
        protParams['value'] = value
        protParams['mode'] = FILTER_LOW_PASS_NO_DECAY
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return FilterVolumesWizard._getListProvider(self, _objs)

    def show(self, form):
        params = self._getParameters(form.protocol)
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            args = {'mode': params['mode'],
                    'highFreq': params['value'],
                    'unit': UNIT_ANGSTROM
                    }

            args['showLowFreq'] = False
            args['showDecay'] = False

            d = BandPassFilterDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar('initialLowPassFilterA', d.samplingRate/d.getHighFreq())

        else:
            dialog.showWarning("Input volumes", "Select volumes first", form.root)
            

class Relion2PartDiameter(RelionPartMaskDiameterWizard):
    _targets = [(ProtRelion2Autopick, ['particleDiameter'])]

    def _getProtocolImages(self, protocol):
        if protocol.useInputReferences():
            return protocol.inputReferences
        else:
            return protocol.inputReferences3D

    def show(self, form, *args):
        prot = form.protocol
        if prot.getInputReferences() is None:
            form.showWarning("Please select the input references first. ")
        else:
            RelionPartMaskDiameterWizard.show(self, form)


class RelionWizMtfSelector(EmWizard):
    """ Simple wizard to select MTF from some of the predefined ones.
    """
    _targets = [(ProtRelionPostprocess, ['mtf']),
                (ProtRelionAssignOpticsGroup, ['mtfFile'])]

    def show(self, form, *args):
        def setPath(fileInfo):
            prot = form.protocol
            varName = 'mtf' if hasattr(prot, 'mtf') else 'mtfFile'
            form.setVar(varName, fileInfo.getPath())
        mtfDir = os.path.join(os.path.dirname(convert.__file__), 'mtfs')
        browser = FileBrowserWindow("Select the one of the predefined MTF files",
                                    form, mtfDir, onSelect=setPath)
        browser.show()


class RelionColorScaleWizard(ColorScaleWizardBase):
    _targets = ColorScaleWizardBase.defineTargets(RelionLocalResViewer)


class RelionWizCtfGroupsDisplay(EmWizard):
    """ Simple wizard distribution of defocus groups.
    """
    _targets = [(ProtRelionClassify2D, ['defocusRange']),
                (ProtRelionClassify3D, ['defocusRange']),
                (ProtRelionRefine3D, ['defocusRange']),
                (ProtRelionInitialModel, ['defocusRange'])]

    def show(self, form, *args):
        prot = form.protocol
        defocusGroups = prot.createDefocusGroups()
        print(defocusGroups)

        plotter = EmPlotter(windowTitle='%d Defocus Groups' % len(defocusGroups),
                            figsize=(8, 6))
        ax = plotter.createSubPlot("", "defocus (A)", "count", 1, 1)

        for group in defocusGroups:
            ax.bar(group.minDefocus, group.count,
                   group.maxDefocus - group.minDefocus,
                   align='edge')

        plotter.show()

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

from collections import OrderedDict

import pwem.emlib.metadata as md

# ----------------- Constants values ------------------------------------------

RELION_HOME = 'RELION_HOME'
RELION_CUDA_LIB = 'RELION_CUDA_LIB'
RELION_PYTHON = 'RELION_PYTHON'

# Supported versions:
V3_1_0 = '3.1.0'
V3_1_1 = '3.1.1'
V3_1_2 = '3.1.2'
V4_0 = '4.0'

MASK_FILL_ZERO = 0
MASK_FILL_NOISE = 1

ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5',
                         '0.2', '0.1', '0.06', '0.03', '0.01', '0.007', '0.004']

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1

# Picking protocols
REF_AVERAGES = 0
REF_VOLUME = REF_BLOBS = 1

# Used in ctf-refinement
FIT_NO = 0
FIT_MICS = 1
FIT_PARTS = 2

# Axis code
AX_X = 0
AX_Y = 1
AX_Z = 2

# Protocol create mask 3d
MASK_AND = 0
MASK_OR = 1
MASK_AND_NOT = 2
MASK_OR_NOT = 3

# Protocol export particles
STACK_NONE = 0
STACK_MULT = 1
STACK_ONE = 2

# Viewer constants
ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CLASSES_ALL = 0
CLASSES_SEL = 1

FSC_CORRECTED = 0
FSC_UNMASKEDMAPS = 1
FSC_MASKEDMAPS = 2
FSC_RANDOMIZED = 3
FSC_ALL = 4

# Color maps
COLOR_JET = 0
COLOR_TERRAIN = 1
COLOR_GIST_EARTH = 2
COLOR_GIST_NCAR = 3
COLOR_GNU_PLOT = 4
COLOR_GNU_PLOT2 = 5
COLOR_OTHER = 6

COLOR_CHOICES = OrderedDict()

COLOR_CHOICES[COLOR_JET] = 'jet'
COLOR_CHOICES[COLOR_TERRAIN] = 'terrain'
COLOR_CHOICES[COLOR_GIST_EARTH] = 'gist_earth'
COLOR_CHOICES[COLOR_GIST_NCAR] = 'gist_ncar'
COLOR_CHOICES[COLOR_GNU_PLOT] = 'gnuplot'
COLOR_CHOICES[COLOR_GNU_PLOT2] = 'gnuplot2'
COLOR_CHOICES[COLOR_OTHER] = 'other'

# This dictionary will be used to map
# between CTFModel properties and Xmipp labels

ACQUISITION_DICT = OrderedDict([
    ("_amplitudeContrast", md.RLN_CTF_Q0),
    ("_sphericalAberration", md.RLN_CTF_CS),
    ("_voltage", md.RLN_CTF_VOLTAGE),
    ("_magnification", md.RLN_CTF_MAGNIFICATION)
])

COOR_DICT = OrderedDict([
    ("_x", md.RLN_IMAGE_COORD_X),
    ("_y", md.RLN_IMAGE_COORD_Y)
])

COOR_EXTRA_LABELS = [
    # Additional autopicking-related metadata
    md.RLN_PARTICLE_AUTOPICK_FOM,
    md.RLN_PARTICLE_CLASS,
    md.RLN_ORIENT_PSI
]

CTF_DICT = OrderedDict([
    ("_defocusU", md.RLN_CTF_DEFOCUSU),
    ("_defocusV", md.RLN_CTF_DEFOCUSV),
    ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
])

CTF_PSD_DICT = OrderedDict([
    ("_psdFile", md.RLN_CTF_IMAGE)
])

CTF_EXTRA_LABELS = [
    md.RLN_CTF_FOM,
    md.RLN_CTF_PHASESHIFT,
    # In Relion the ctf also contains acquisition information
    md.RLN_CTF_Q0,
    md.RLN_CTF_CS,
    md.RLN_CTF_VOLTAGE,
    md.RLN_CTF_MAGNIFICATION,
    md.RLN_CTF_DETECTOR_PIXEL_SIZE
]

# Some extra labels
IMAGE_EXTRA_LABELS = [
    md.RLN_SELECT_PARTICLES_ZSCORE,
    md.RLN_IMAGE_FRAME_NR
]

ALIGNMENT_DICT = OrderedDict([
    ("_rlnOriginX", md.RLN_ORIENT_ORIGIN_X),
    ("_rlnOriginY", md.RLN_ORIENT_ORIGIN_Y),
    ("_rlnOriginZ", md.RLN_ORIENT_ORIGIN_Z),
    ("_rlnAngleRot", md.RLN_ORIENT_ROT),
    ("_rlnAngleTilt", md.RLN_ORIENT_TILT),
    ("_rlnAnglePsi", md.RLN_ORIENT_PSI)
])


PARTICLE_EXTRA_LABELS = [
    'rlnCtfBfactor',
    'rlnCtfScalefactor',
    'rlnPhaseShift',
    'rlnGroupName',
    'rlnGroupNumber',
    'rlnRandomSubset'
]

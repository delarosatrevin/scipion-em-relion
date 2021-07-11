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

import pyworkflow.utils as pwutils
import pwem

from .constants import *


__version__ = '4.0'
_logo = "relion_logo.jpg"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016', 'Zivanov2018']


class Plugin(pwem.Plugin):
    _homeVar = RELION_HOME
    _supportedVersions = [V3_1_0, V3_1_1, V3_1_2, V4_0]
    _url = "https://github.com/scipion-em/scipion-em-relion"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(RELION_HOME, 'relion-%s' % V4_0)
        cls._defineVar(RELION_CUDA_LIB, pwem.Config.CUDA_LIB)
        cls._defineVar(RELION_PYTHON, None)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Relion. """
        environ = pwutils.Environ(os.environ)
        binPath = os.pathsep.join([cls.getHome('bin'),
                                   pwem.Config.MPI_BINDIR])
        libPath = os.pathsep.join([cls.getHome('lib'),
                                   cls.getHome('lib64'),
                                   pwem.Config.MPI_LIBDIR])

        if binPath not in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath
                            }, position=pwutils.Environ.BEGIN)

        # Get Relion CUDA library path if defined
        cudaLib = cls.getVar(RELION_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)

        if 'RELION_MPI_LIB' in os.environ:
            environ.addLibrary(os.environ['RELION_MPI_LIB'])

        if 'RELION_MPI_BIN' in os.environ:
            environ.set('PATH', os.environ['RELION_MPI_BIN'],
                        position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def IS_GT31(cls):
        return cls.getActiveVersion().startswith('4')

    @classmethod
    def defineBinaries(cls, env):
        relion_commands = [('cmake -DGUI=OFF -DCMAKE_INSTALL_PREFIX=./ .', []),
                           ('make -j %d' % env.getProcessors(),
                            ['bin/relion_refine'])]

        for v in cls._supportedVersions:
            env.addPackage('relion', version=v,
                           url='https://github.com/3dem/relion/archive/%s.tar.gz' % v,
                           commands=relion_commands,
                           updateCuda=True,
                           default=v == V4_0)

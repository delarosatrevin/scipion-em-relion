# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
import pyworkflow.em
import pyworkflow.utils as pwutils

from .constants import RELION_HOME, V2_0, V2_1


_logo = "relion_logo.png"


class Plugin(pyworkflow.em.Plugin):
    _homeVar = RELION_HOME
    _pathVars = [RELION_HOME]
    _supportedVersions = [V2_0, V2_1]

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Relion. """

        environ = pwutils.Environ(os.environ)
        binPath = cls.getHome('bin')
        libPath = cls.getHome('lib') + ":" + cls.getHome('lib64')

        if not binPath in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath,
                            'SCIPION_MPI_FLAGS': os.environ.get('RELION_MPI_FLAGS', ''),
                            }, position=pwutils.Environ.BEGIN)

        # Take Scipion CUDA library path
        cudaLib = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def isVersion2Active(cls):
        return cls.getActiveVersion().startswith("2.")


pyworkflow.em.Domain.registerPlugin(__name__)
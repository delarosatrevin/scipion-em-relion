# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
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
"""
New conversion functions dealing with Relion3.1 new star files format.
"""
import os
import io
import numpy as np
from collections import OrderedDict
from emtable import Table


from pwem.constants import ALIGN_NONE, ALIGN_PROJ, ALIGN_2D, ALIGN_3D
from pwem.objects import (Micrograph, SetOfMicrographsBase, SetOfMovies,
                          Particle, CTFModel, Acquisition, Transform, Coordinate)
import pwem.convert.transformations as tfs

from .convert_base import WriterBase, ReaderBase
from .convert_utils import (convertBinaryFiles, locationToRelion,
                            relionToLocation)
from relion.constants import PARTICLE_EXTRA_LABELS


def getPixelSizeLabel(imageSet):
    """ Return the proper label for pixel size. """
    if (isinstance(imageSet, SetOfMicrographsBase)
            or isinstance(imageSet, Micrograph)):
        return 'rlnMicrographPixelSize'
    else:
        return 'rlnImagePixelSize'


class OpticsGroups:
    """ Store information about optics groups in an indexable way.
    Existing groups can be accessed by number of name.
    """

    def __init__(self, opticsTable):
        self.__fromTable(opticsTable)

    def __fromTable(self, opticsTable):
        self._dict = OrderedDict()
        # Also allow indexing by name
        self._dictName = OrderedDict()
        # Map optics rows both by name and by number
        for og in opticsTable:
            self.__store(og)

    def __store(self, og):
        self._dict[og.rlnOpticsGroup] = og
        self._dictName[og.rlnOpticsGroupName] = og

    def __getitem__(self, item):
        if isinstance(item, int):
            return self._dict[item]
        elif isinstance(item, str):
            return self._dictName[item]
        raise TypeError("Unsupported type '%s' of item '%s'"
                        % (type(item), item))

    def __contains__(self, item):
        return item in self._dict or item in self._dictName

    def __iter__(self):
        """ Iterate over all optics groups. """
        return iter(self._dict.values())

    def __len__(self):
        return len(self._dict)

    def __str__(self):
        return self.toString()

    def first(self):
        """ Return first optics group. """
        return next(iter(self._dict.values()))

    def update(self, ogId, **kwargs):
        og = self.__getitem__(ogId)
        newOg = og._replace(**kwargs)
        self.__store(newOg)
        return newOg

    def updateAll(self, **kwargs):
        """ Update all Optics Groups with these values. """
        missing = {k: v for k, v in kwargs.items() if not self.hasColumn(k)}
        existing = {k: v for k, v in kwargs.items() if self.hasColumn(k)}

        self.addColumns(**missing)

        for og in self:
            self.update(og.rlnOpticsGroup, **existing)

    def add(self, newOg):
        self.__store(newOg)

    def hasColumn(self, colName):
        return hasattr(self.first(), colName)

    def addColumns(self, **kwargs):
        """ Add new columns with default values (type inferred from it). """
        items = self.first()._asdict().items()
        cols = [Table.Column(k, type(v)) for k, v in items]

        for k, v in kwargs.items():
            cols.append(Table.Column(k, type(v)))

        t = Table(columns=cols)

        for og in self._dict.values():
            values = og._asdict()
            values.update(kwargs)
            t.addRow(**values)

        self.__fromTable(t)

    @staticmethod
    def fromStar(starFilePath):
        """ Create an OpticsGroups from a given STAR file.
        """
        return OpticsGroups(Table(fileName=starFilePath, tableName='optics'))

    @staticmethod
    def fromString(stringValue):
        """ Create an OpticsGroups from string content (STAR format)
        """
        f = io.StringIO(stringValue)
        t = Table()
        t.readStar(f, tableName='optics')
        return OpticsGroups(t)

    @staticmethod
    def fromImages(imageSet):
        acq = imageSet.getAcquisition()
        params = {'rlnImageSize': imageSet.getXDim(),
                  getPixelSizeLabel(imageSet): imageSet.getSamplingRate()}
        if isinstance(imageSet, SetOfMovies):
            params['rlnMicrographOriginalPixelSize'] = imageSet.getSamplingRate()
        try:
            og = OpticsGroups.fromString(acq.opticsGroupInfo.get())
            # always update sampling and image size from the set
            og.updateAll(**params)
            return og
        except:
            params.update({
                'rlnVoltage': acq.getVoltage(),
                'rlnSphericalAberration': acq.getSphericalAberration(),
                'rlnAmplitudeContrast': acq.getAmplitudeContrast(),
            })
            return OpticsGroups.create(**params)

    @staticmethod
    def create(**kwargs):
        opticsString1 = """

# version 40001

data_optics

loop_ 
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnMicrographOriginalPixelSize #3
_rlnVoltage #4
_rlnSphericalAberration #5
_rlnAmplitudeContrast #6
_rlnImageSize #7
_rlnImageDimensionality #8
opticsGroup1            1      1.000000   300.000000     2.700000     0.100000     256            2
        """

        og = OpticsGroups.fromString(opticsString1)
        fog = og.first()
        newColumns = {k: v for k, v in kwargs.items() if not hasattr(fog, k)}
        og.addColumns(**newColumns)
        og.update(1, **kwargs)
        return og

    def _write(self, f):
        # Create columns from the first row
        items = self.first()._asdict().items()
        cols = [Table.Column(k, type(v)) for k, v in items]
        t = Table(columns=cols)
        for og in self._dict.values():
            t.addRow(*og)
        t.writeStar(f, tableName='optics')

    def toString(self):
        """ Return a string (STAR format) with the current optics groups.
        """
        f = io.StringIO()
        self._write(f)
        result = f.getvalue()
        f.close()

        return result

    def toStar(self, starFile):
        """ Write current optics groups to a given file.
        """
        self._write(starFile)

    def toImages(self, imageSet):
        """ Store the optics groups information in the image acquisition.
        """
        imageSet.getAcquisition().opticsGroupInfo.set(self.toString())


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfImages subclasses
    into Relion>3.1 star files (and binaries if conversion needed).
    """

    def writeSetOfMovies(self, moviesIterable, starFile, **kwargs):
        self._writeSetOfMoviesOrMics(moviesIterable, starFile,
                                     'movies', 'rlnMicrographMovieName',
                                     **kwargs)

    def writeSetOfMicrographs(self, micsIterable, starFile, **kwargs):
        self._writeSetOfMoviesOrMics(micsIterable, starFile,
                                     'micrographs', 'rlnMicrographName',
                                     **kwargs)

    def _writeSetOfMoviesOrMics(self, imgIterable,
                                starFile, tableName, imgLabelName, **kwargs):
        """ This function can be used to write either movies or micrographs
        star files. Input can be any iterable of these type of images (e.g
        set, list, etc).
        """
        # Process the first item and create the table based
        # on the generated columns
        self._imgLabelName = imgLabelName
        self._postprocessImageRow = kwargs.get('postprocessImageRow', None)
        self._prefix = tableName[:3]

        micRow = OrderedDict()
        micRow[imgLabelName] = ''  # Just to add label, proper value later
        iterMics = iter(imgIterable)
        mic = next(iterMics)
        if self._optics is None:
            self._optics = OpticsGroups.fromImages(mic)
        self._imageSize = mic.getXDim()
        self._setCtf = mic.hasCTF()

        extraLabels = kwargs.get('extraLabels', [])
        self._extraLabels = [l for l in extraLabels
                             if mic.hasAttribute('_%s' % l)]

        self._micToRow(mic, micRow)
        if self._postprocessImageRow:
            self._postprocessImageRow(mic, micRow)

        micsTable = self._createTableFromDict(micRow)

        while mic is not None:
            micRow[imgLabelName] = self._convert(mic)
            self._micToRow(mic, micRow)

            if self._postprocessImageRow:
                self._postprocessImageRow(mic, micRow)

            micsTable.addRow(**micRow)
            mic = next(iterMics, None)

        with open(starFile, 'w') as f:
            f.write("# Star file generated with Scipion\n")
            f.write("# version 40001\n")
            self._optics.toStar(f)
            f.write("# version 40001\n")
            micsTable.writeStar(f, tableName=tableName)

    def _objToRow(self, obj, row, attributes):
        """ Set some attributes from the object to the row.
        For performance reasons, it is not validated that each attribute
        is already in the object, so it should be validated before.
        """
        for attr in attributes:
            row[attr] = obj.getAttributeValue('_%s' % attr)

    def _micToRow(self, mic, row):
        WriterBase._micToRow(self, mic, row)

        # Set CTF values
        if self._setCtf:
            self._ctfToRow(mic.getCTF(), row)

        # Set additional labels if present
        self._objToRow(mic, row, self._extraLabels)
        row['rlnOpticsGroup'] = mic.getAttributeValue('_rlnOpticsGroup', 1)

    def _align2DToRow(self, alignment, row):
        matrix = alignment.getMatrix()
        shifts = tfs.translation_from_matrix(matrix)
        shifts *= self._pixelSize
        angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
        row['rlnOriginXAngst'], row['rlnOriginYAngst'] = shifts[:2]
        row['rlnAnglePsi'] = -(angles[0] + angles[2])

    def _alignProjToRow(self, alignment, row):
        matrix = np.linalg.inv(alignment.getMatrix())
        shifts = -tfs.translation_from_matrix(matrix)
        shifts *= self._pixelSize
        angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
        row['rlnOriginXAngst'], row['rlnOriginYAngst'], row['rlnOriginZAngst'] = shifts
        row['rlnAngleRot'], row['rlnAngleTilt'], row['rlnAnglePsi'] = angles

    def _partToRow(self, part, row):
        row['rlnImageId'] = part.getObjId()

        # Add coordinate information
        coord = part.getCoordinate()
        if coord is not None:
            x, y = coord.getPosition()
            row['rlnCoordinateX'] = x
            row['rlnCoordinateY'] = y
            # Add some specific coordinate attributes
            self._objToRow(coord, row, self._coordLabels)
            micName = coord.getMicName()
            if micName:
                row['rlnMicrographName'] = str(micName.replace(" ", ""))
            else:
                if coord.getMicId():
                    row['rlnMicrographName'] = str(coord.getMicId())

        index, fn = part.getLocation()
        if self.outputStack:
            row['rlnOriginalParticleName'] = locationToRelion(index, fn)
            index, fn = self._counter, self._relOutputStack
            if self._counter > 0:
                self._ih.convert(part, (index, self.outputStack))
        else:
            if self.outputDir is not None:
                fn = self._filesDict.get(fn, fn)

        row['rlnImageName'] = locationToRelion(index, fn)

        # Set CTF values
        if self._setCtf:
            self._ctfToRow(part.getCTF(), row)

        # Set alignment if necessary
        if self._setAlign:
            self._setAlign(part.getTransform(), row)

        # Set additional labels if present
        self._objToRow(part, row, self._extraLabels)

        # Add now the new Optics Group stuff
        row['rlnOpticsGroup'] = part.getAttributeValue('_rlnOpticsGroup', 1)

        self._counter += 1

    def writeSetOfParticles(self, partsSet, starFile, **kwargs):
        # Process the first item and create the table based
        # on the generated columns
        self.update(['rootDir', 'outputDir', 'outputStack'], **kwargs)

        self._optics = OpticsGroups.fromImages(partsSet)
        partRow = OrderedDict()
        firstPart = partsSet.getFirstItem()

        # Convert binaries if required
        if self.outputStack:
            self._relOutputStack = os.path.relpath(self.outputStack,
                                                   os.path.dirname(starFile))
        if self.outputDir is not None:
            forceConvert = kwargs.get('forceConvert', False)
            self._filesDict = convertBinaryFiles(partsSet, self.outputDir,
                                                 forceConvert=forceConvert)

        # Compute some flags from the first particle...
        # when flags are True, some operations will be applied to all particles
        self._preprocessImageRow = kwargs.get('preprocessImageRow', None)

        self._setCtf = kwargs.get('writeCtf', True) and firstPart.hasCTF()

        alignType = kwargs.get('alignType', partsSet.getAlignment())

        if alignType == ALIGN_2D:
            self._setAlign = self._align2DToRow
        elif alignType == ALIGN_PROJ:
            self._setAlign = self._alignProjToRow
        elif alignType == ALIGN_3D:
            raise NotImplementedError(
                "3D alignment conversion for Relion not implemented. "
                "It seems the particles were generated with an incorrect "
                "alignment type. You may either re-launch the protocol that "
                "generates the particles with angles or set 'Consider previous"
                " alignment?' to No")
        elif alignType == ALIGN_NONE:
            self._setAlign = None
        else:
            raise TypeError("Invalid value for alignType: %s" % alignType)

        extraLabels = kwargs.get('extraLabels', [])
        extraLabels.extend(PARTICLE_EXTRA_LABELS)
        self._extraLabels = [l for l in extraLabels
                             if firstPart.hasAttribute('_%s' % l)]

        coord = firstPart.getCoordinate()
        self._coordLabels = []
        if coord is not None:
            self._coordLabels = [l for l in ['rlnClassNumber',
                                             'rlnAutopickFigureOfMerit',
                                             'rlnAnglePsi']
                                 if coord.hasAttribute('_%s' % l)]

        self._postprocessImageRow = kwargs.get('postprocessImageRow', None)

        self._imageSize = firstPart.getXDim()
        self._pixelSize = firstPart.getSamplingRate() or 1.0

        self._counter = 0  # Mark first conversion as special one
        firstPart.setAcquisition(partsSet.getAcquisition())
        self._partToRow(firstPart, partRow)

        if self._postprocessImageRow:
            self._postprocessImageRow(firstPart, partRow)

        partsTable = self._createTableFromDict(partRow)
        partsTable.addRow(**partRow)

        with open(starFile, 'w') as f:
            # Write particles table
            f.write("# Star file generated with Scipion\n")
            f.write("\n# version 40001\n")
            self._optics.toStar(f)
            f.write("# version 40001\n")
            # Write header first
            partsWriter = Table.Writer(f)
            partsWriter.writeTableName('particles')
            partsWriter.writeHeader(partsTable.getColumns())
            # Write all rows
            for part in partsSet:
                self._partToRow(part, partRow)
                if self._postprocessImageRow:
                    self._postprocessImageRow(part, partRow)
                partsWriter.writeRowValues(partRow.values())
                # partsTable.writeStarLine(f, partRow.values())


class Reader(ReaderBase):
    ALIGNMENT_LABELS = [
        "rlnOriginXAngst",
        "rlnOriginYAngst",
        "rlnOriginZAngst",
        "rlnAngleRot",
        "rlnAngleTilt",
        "rlnAnglePsi",
    ]

    CTF_LABELS = [
        "rlnDefocusU",
        "rlnDefocusV",
        "rlnDefocusAngle",
        "rlnCtfAstigmatism",
        "rlnCtfFigureOfMerit",
        "rlnCtfMaxResolution"
    ]

    COORD_LABELS = [
        "rlnCoordinateX",
        "rlnCoordinateY",
        "rlnMicrographName",
        # extra labels below
        "rlnAutopickFigureOfMerit",
        "rlnClassNumber",
        "rlnAnglePsi"
    ]

    def __init__(self, **kwargs):
        """
        """
        ReaderBase.__init__(self, **kwargs)

    def readSetOfParticles(self, starFile, partSet, **kwargs):
        """ Convert a star file into a set of particles.

        Params:
            starFile: the filename of the star file
            partsSet: output particles set

        Keyword Arguments:
            blockName: The name of the data block (default particles)
            alignType: alignment type
            removeDisabled: Remove disabled items

        """
        self._preprocessImageRow = kwargs.get('preprocessImageRow', None)
        self._alignType = kwargs.get('alignType', ALIGN_NONE)

        self._postprocessImageRow = kwargs.get('postprocessImageRow', None)

        self._optics = OpticsGroups.fromStar(starFile)

        self._pixelSize = getattr(self._optics.first(),
                                  'rlnImagePixelSize', 1.0)
        self._invPixelSize = 1. / self._pixelSize

        partsReader = Table.Reader(starFile, tableName='particles')

        firstRow = partsReader.getRow()
        self._setClassId = hasattr(firstRow, 'rlnClassNumber')
        self._setCtf = partsReader.hasAllColumns(self.CTF_LABELS[:3])
        self._setCoord = partsReader.hasAllColumns(self.COORD_LABELS[:3])
        particle = Particle()

        if self._setCtf:
            particle.setCTF(CTFModel())

        self._setAcq = kwargs.get("readAcquisition", True)
        if self._setAcq:
            acq = Acquisition()
            self.rowToAcquisition(self._optics.first(), acq)
            acq.setMagnification(kwargs.get('magnification', 10000))
            partSet.setAcquisition(acq)
        else:
            # readAcquisition=False ONLY during import particles
            # overwrite pixel size and optics
            self._pixelSize = kwargs.get('pixelSize', self._pixelSize)
            acq = partSet.getAcquisition()
            self._optics.updateAll(
                rlnVoltage=acq.getVoltage(),
                rlnSphericalAberration=acq.getSphericalAberration(),
                rlnAmplitudeContrast=acq.getAmplitudeContrast(),
                rlnImagePixelSize=self._pixelSize
            )

        extraLabels = kwargs.get('extraLabels', []) + PARTICLE_EXTRA_LABELS
        self.createExtraLabels(particle, firstRow, extraLabels)

        self._rowToPart(firstRow, particle)
        partSet.setSamplingRate(self._pixelSize)
        self._optics.toImages(partSet)
        partSet.append(particle)

        for row in partsReader:
            self._rowToPart(row, particle)
            partSet.append(particle)

        partSet.setHasCTF(self._setCtf)
        partSet.setAlignment(self._alignType)

    def _rowToPart(self, row, particle):
        particle.setObjId(getattr(row, 'rlnImageId', None))

        if self._preprocessImageRow:
            self._preprocessImageRow(particle, row)

        # Decompose Relion filename
        index, filename = relionToLocation(row.rlnImageName)
        particle.setLocation(index, filename)

        if self._setClassId:
            particle.setClassId(row.rlnClassNumber)

        if self._setCtf:
            self.rowToCtf(row, particle.getCTF())

        self.setParticleTransform(particle, row)
        self.setExtraLabels(particle, row)

        # TODO: coord extra labels, partId, micId,
        if self._setCoord:
            coord = self.rowToCoord(row)
            particle.setCoordinate(coord)

        if self._postprocessImageRow:
            self._postprocessImageRow(particle, row)

    @staticmethod
    def rowToCoord(row):
        """ Create a Coordinate from the row. """
        coord = Coordinate()
        coord.setPosition(row.rlnCoordinateX,
                          row.rlnCoordinateY)
        coord.setMicName(row.rlnMicrographName)

        return coord

    @staticmethod
    def rowToCtf(row, ctf):
        """ Create a CTFModel from the row. """
        ctf.setDefocusU(row.rlnDefocusU)
        ctf.setDefocusV(row.rlnDefocusV)
        ctf.setDefocusAngle(row.rlnDefocusAngle)
        ctf.setResolution(row.get('rlnCtfMaxResolution', 0))
        ctf.setFitQuality(row.get('rlnCtfFigureOfMerit', 0))

        if hasattr(row, 'rlnPhaseShift'):
            ctf.setPhaseShift(row.rlnPhaseShift)
        ctf.standardize()

        if hasattr(row, 'rlnCtfImage'):
            ctf.setPsdFile(row.rlnCtfImage)

    @staticmethod
    def rowToAcquisition(optics, acq):
        acq.setAmplitudeContrast(optics.rlnAmplitudeContrast)
        acq.setSphericalAberration(optics.rlnSphericalAberration)
        acq.setVoltage(optics.rlnVoltage)

    def setParticleTransform(self, particle, row):
        """ Set the transform values from the row. """

        if ((self._alignType == ALIGN_NONE) or
                not row.hasAnyColumn(self.ALIGNMENT_LABELS)):
            self.setParticleTransform = self.__setParticleTransformNone
        else:
            # Ensure the Transform object exists
            self._angles = np.zeros(3)
            self._shifts = np.zeros(3)

            particle.setTransform(Transform())

            if self._alignType == ALIGN_2D:
                self.setParticleTransform = self.__setParticleTransform2D
            elif self._alignType == ALIGN_PROJ:
                self.setParticleTransform = self.__setParticleTransformProj
            else:
                raise TypeError("Unexpected alignment type: %s"
                                % self._alignType)

        # Call again the modified function
        self.setParticleTransform(particle, row)

    def __setParticleTransformNone(self, particle, row):
        particle.setTransform(None)

    def __setParticleTransform2D(self, particle, row):
        angles = self._angles
        shifts = self._shifts
        ips = self._invPixelSize

        def _get(label):
            return float(getattr(row, label, 0.))

        shifts[0] = _get('rlnOriginXAngst') * ips
        shifts[1] = _get('rlnOriginYAngst') * ips
        angles[2] = -_get('rlnAnglePsi')
        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = shifts[:3]
        particle.getTransform().setMatrix(M)

    def __setParticleTransformProj(self, particle, row):
        angles = self._angles
        shifts = self._shifts
        ips = self._invPixelSize

        def _get(label):
            return float(getattr(row, label, 0.))

        shifts[0] = _get('rlnOriginXAngst') * ips
        shifts[1] = _get('rlnOriginYAngst') * ips
        shifts[2] = _get('rlnOriginZAngst') * ips

        angles[0] = _get('rlnAngleRot')
        angles[1] = _get('rlnAngleTilt')
        angles[2] = _get('rlnAnglePsi')

        radAngles = -np.deg2rad(angles)

        # TODO: jmrt: Maybe we should test performance and consider if keeping
        # TODO: the matrix and not creating one everytime will make things faster
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
        particle.getTransform().setMatrix(M)

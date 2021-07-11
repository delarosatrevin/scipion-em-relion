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
Utility functions for conversions that will be used from both
newer Relion3.1 routines and old ones.
"""

import os
import math
from emtable import Table

import pyworkflow.utils as pwutils
from pwem.constants import NO_INDEX
from pwem.emlib.image import ImageHandler

from relion import Plugin


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    if index != NO_INDEX:
        return "%06d@%s" % (index, filename)

    return filename


def relionToLocation(filename):
    """ Return a location (index, filename) given
    a Relion filename with the index@filename structure. """
    if '@' in filename:
        indexStr, fn = filename.split('@')
        return int(indexStr), str(fn)
    else:
        return NO_INDEX, str(filename)


def convertBinaryFiles(imgSet, outputDir, extension='mrcs', forceConvert=False):
    """ Convert binary images files to a format read by Relion.
    Or create links if there is no need to convert the binary files.

    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
        extension: extension accepted by the program
        forceConvert: if True, the files will be converted and no root will be used
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesDict = {}
    ih = ImageHandler()
    outputRoot = outputDir if forceConvert else os.path.join(outputDir, 'input')
    # Get the extension without the dot
    stackFiles = imgSet.getFiles()
    ext = pwutils.getExt(next(iter(stackFiles)))[1:]
    rootDir = pwutils.commonPath(list(stackFiles))

    def getUniqueFileName(fn, extension):
        """ Get an unique file for either link or convert files.
        It is possible that the base name overlap if they come
        from different runs. (like particles.mrcs after relion preprocess)
        """
        newFn = os.path.join(outputRoot, pwutils.replaceBaseExt(fn, extension))
        newRoot = pwutils.removeExt(newFn)

        values = filesDict.values()
        counter = 1

        while newFn in values:
            counter += 1
            newFn = '%s_%05d.%s' % (newRoot, counter, extension)

        return newFn

    def createBinaryLink(fn):
        """ Just create a link named .mrcs to Relion understand
        that it is a binary stack file and not a volume.
        """
        newFn = getUniqueFileName(fn, extension)
        if not os.path.exists(newFn):
            pwutils.createLink(fn, newFn)
            print("   %s -> %s" % (newFn, fn))
        return newFn

    def convertStack(fn):
        """ Convert from a format that is not read by Relion
        to an spider stack.
        """
        newFn = getUniqueFileName(fn, 'mrcs')
        ih.convertStack(fn, newFn)
        print("   %s -> %s" % (newFn, fn))
        return newFn

    def replaceRoot(fn):
        """ Link create to the root folder, so just replace that
        in the name, no need to do anything else.
        """
        return fn.replace(rootDir, outputRoot)

    if forceConvert:
        print("convertBinaryFiles: forceConvert = True")
        mapFunc = convertStack
    elif ext == extension:
        print("convertBinaryFiles: creating soft links.")
        print("   Root: %s -> %s" % (outputRoot, rootDir))
        mapFunc = replaceRoot
        pwutils.createAbsLink(os.path.abspath(rootDir), outputRoot)
    elif ext == 'mrc' and extension == 'mrcs':
        print("convertBinaryFiles: creating soft links (mrcs -> mrc).")
        mapFunc = createBinaryLink
    elif ext.endswith('hdf'):  # assume eman .hdf format
        print("convertBinaryFiles: converting stacks. (%s -> %s)"
              % (extension, ext))
        mapFunc = convertStack
    else:
        mapFunc = None

    if mapFunc is not None:
        pwutils.makePath(outputRoot)
        for fn in stackFiles:
            newFn = mapFunc(fn)  # convert or link
            filesDict[fn] = newFn  # map new filename

    return filesDict


def convertBinaryVol(vol, outputDir):
    """ Convert binary volume to a format read by Relion.
    Params:
        vol: input volume object to be converted.
        outputDir: where to put the converted file(s)
    Return:
        new file name of the volume (converted or not).
    """

    ih = ImageHandler()
    fn = vol.getFileName()

    if not fn.endswith('.mrc'):
        newFn = os.path.join(outputDir, pwutils.replaceBaseExt(fn, 'mrc'))
        ih.convert(fn, newFn)
        return newFn

    return fn


def convertMask(img, outputPath, newPix=None,
                newDim=None, threshold=True, invert=False):
    """ Convert mask to mrc format read by Relion.
    Params:
        img: input image to be converted.
        outputPath: it can be either a directory or a file path.
            If it is a directory, the output name will be inferred from input
            and put into that directory. If it is not a directory,
            it is assumed is the output filename.
        newPix: output pixel size (equals input if None)
        newDim: output box size
    Return:
        new file name of the mask.
    """
    index, filename = img.getLocation()
    imgFn = locationToRelion(index, filename)
    inPix = img.getSamplingRate()

    if os.path.isdir(outputPath):
        outFn = os.path.join(outputPath, pwutils.replaceBaseExt(imgFn, 'mrc'))
    else:
        outFn = outputPath

    if invert:
        # unfortunately relion does not allow to use
        # add_constant and multiply_constant together
        # so let us multiply first
        params = '--i %s --o %s --angpix %0.5f' % (imgFn, outFn, inPix)
        params += ' --multiply_constant -1'
        pwutils.runJob(None, 'relion_image_handler', params, env=Plugin.getEnviron())
        imgFn = outFn

    params = '--i %s --o %s --angpix %0.5f' % (imgFn, outFn, inPix)

    if newPix is not None and not math.isclose(newPix, inPix, abs_tol=0.001):
        # be careful with this rescale param because it may
        # introduce some artefacts if the sampling is the same as the input
        params += ' --rescale_angpix %0.5f' % newPix

    if newDim is not None:
        params += ' --new_box %d' % newDim

    if threshold:
        params += ' --threshold_above 1 --threshold_below 0'

    if invert:
        # unfortunately relion does not allow to use
        # add_constant and multiply_constant together
        params += ' --add_constant 1'

    pwutils.runJob(None, 'relion_image_handler', params, env=Plugin.getEnviron())
    return outFn


def relativeFromFileName(imgRow, prefixPath):
    """ Remove some prefix from filename in row. """
    index, imgPath = relionToLocation(imgRow['rlnImageName'])
    newImgPath = os.path.relpath(imgPath, prefixPath)
    imgRow['rlnImageName'] = locationToRelion(index, newImgPath)


def getVolumesFromPostprocess(postStar):
    """ Return the filenames of half1, half2 and mask from
    a given postprocess.star file.
    """
    table = Table(fileName=postStar, tableName='general')
    row = table[0]
    return (row.rlnUnfilteredMapHalf1,
            row.rlnUnfilteredMapHalf2,
            row.rlnMaskName)

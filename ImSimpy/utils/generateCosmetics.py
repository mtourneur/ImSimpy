"""
CCD Cosmetics Maps
==================

A simple script to generate a random location cosmetics map.

:requires: NumPy

:author: Sami-Matias Niemi
:contact: s.niemi@ucl.ac.uk
"""
import os, datetime
import numpy as np
from astropy.io import fits

__version__ = 1.0


def HotPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/HotPixs_vis.fits', Type='random',hotpixels=500, hotmin=70000,hotmax=300000,xsize=4096, ysize=4096):
    """
    Generates a map of hot pixels with random locations.

    :return: None
    """

    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if Type == 'random':

        hotpix_map=np.zeros([xsize,ysize])
        xcoords = np.random.random((hotpixels))*xsize
        ycoords = np.random.random((hotpixels))*ysize
        hotvalues = np.random.randint(hotmin, hotmax, hotpixels)

        for x, y, hot in zip(xcoords, ycoords,hotvalues):
            hotpix_map[int(x),int(y)]=hot

    #save to file
    hdu=fits.PrimaryHDU()

    #Normalised surface
    hdu.data=hotpix_map.astype(np.float32)

    if Type == 'random': hdu.header.set('HOTPIXS', 'Map of %d Hot pixels generated randomly' % int(hotpixels) )

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)


def DeadPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/DeadPixs_vis.fits', Type='random',deadpixels=2000, xsize=4096, ysize=4096):
    """
    Generates a map of hot pixels with random locations.

    :return: None
    """

    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if Type == 'random':
        deadpix_map=np.ones([xsize,ysize])
        xcoords = np.random.random((deadpixels))*xsize
        ycoords = np.random.random((deadpixels))*ysize

        for x, y in zip(xcoords, ycoords):
            deadpix_map[int(x),int(y)]=0


    #save to file
    hdu=fits.PrimaryHDU()

    #Normalised surface
    hdu.data=deadpix_map.astype(np.float32)

    if Type == 'random': hdu.header.set('DEADPIX', 'Map of %d dead pixels generated randomly' % int(deadpixels) )

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)


def makeCosmetics(output_dir='Calibration/COLIBRI/',suffix=['_vis','_nir'],xsize=[4096,2048],ysize=[4096,2048],hotpixels=[500,500],hotmin=[70000,5000],hotmax=[300000,40000],deadpixels=[2000,500]):
    # Make hot pixels and dead pixels maps
    from ImSimpy.utils.generateCosmetics import HotPixs, DeadPixs

    for i in range(len(suffix)):
        HotPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/%s/HotPixs%s.fits' % (output_dir,suffix[i]), xsize=xsize[i],ysize=ysize[i],Type='random',hotpixels=hotpixels[i],hotmin=hotmin[i],hotmax=hotmax[i])

        DeadPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/%s/DeadPixs%s.fits' % (output_dir,suffix[i]), xsize=xsize[i],ysize=ysize[i],Type='random',deadpixels=deadpixels[i])
    
    
if __name__ == '__main__':

    #Visible
    HotPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/HotPixs_vis.fits', Type='random',hotpixels=500, hotmin=70000,hotmax=300000,xsize=4096, ysize=4096)
    DeadPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/DeadPixs_vis.fits', Type='random',deadpixels=2000, xsize=4096, ysize=4096)

    #NIR
    HotPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/HotPixs_nir.fits', Type='random',hotpixels=500, hotmin=5000,hotmax=40000,xsize=2048, ysize=2048)
    DeadPixs(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Cosmetics/DeadPixs_nir.fits', Type='random',deadpixels=500, xsize=2048, ysize=2048)

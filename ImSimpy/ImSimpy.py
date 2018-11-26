# -*- coding: utf-8 -*-

"""Main module."""

"""
Optical/NIR Telescope Image Simulator
=============================================

This file contains an image simulatorfor optical/NIR telescope

Contact Information
-------------------
mail: david.corre@lam.fr

"""
import os, sys, datetime, math
from collections import OrderedDict
import hjson
import scipy
from scipy.ndimage import interpolation,zoom
from scipy import ndimage
from astropy.io import fits
import numpy as np
import numexpr as ne
from ImSimpy.utils.generateCalib import DarkCurrent
from ImSimpy.utils import PSFUtils
from scipy.signal import fftconvolve

convolution = fftconvolve

class ImageSimulator():
    """
    Image Simulator for optical/NIR telescope

    """

    def __init__(self, path=os.getenv('ImSimpy_DIR')+'\\ImSimpy',configFile=os.getenv('ImSimpy_DIR')+'/ImSimpy/configFiles/default_input.hjson',name_telescope='default',seed=None, debug=False,random=False):
        """
        Class Constructor.

        :configFile : name of the cofig file 
        :seed: value used to set the random number generator (default: None --> use current time).
        """
        self.path = path
        self.configfile = configFile
        self.name_telescope = name_telescope

        if seed != None:
            seed = int(seed)
            np.random.seed(seed=seed)  #fix the seed
            print ('Setting the random number generator seed to {}'.format(seed))
        else:
            np.random.seed()
            print ('Setting the random number generator seed: current time')
        self.debug=bool(debug)
        self.random=bool(random)

        #load instrument model, these values are also stored in the FITS header
        self.information = OrderedDict()

        #update settings with defaults. It will be used if they are missing if the config file
        self.information.update(dict(psfoversampling=1.0,
                                     xsize=4096,
                                     ysize=4096,
                                     FWC=350000,
                                     DC=0.001,
                                     RN=8.0,
                                     bias=500.0,
                                     gain=3.0,
                                     zeropoint=23.5,
                                     Nexp=1,
                                     exptime=565.0,
                                     readouttime=5.,
                                     sky_brightness=19.3,
                                     RA=123.0,
                                     DEC=45.0,
                                     mode='same'))


    def readConfigs(self):
        """
        Reads the config file information using configParser.
        """
        # Load the input file in hjson format into a dictionary
        with open(self.configfile,encoding='utf-8') as f:
            self.config=hjson.load(f)

    def etc(self,config_type):
        """ Execute the Exposure Time Calculator to get some information (zeropoint, grb mag,...) """

        try:
           from pyETC.pyETC import etc
        except ValueError:
           print ('Package ETC not found, you have to install it')
        if config_type =='file': etc_info=etc(configFile=self.configfile, config_type=config_type,name_telescope=self.name_telescope)
        elif config_type =='data': etc_info=etc(configFile=self.config, config_type=config_type,name_telescope=self.name_telescope)
        etc_info.sim()
        #Update the config file 
        self.config['camera_type']=etc_info.information['cameras'][etc_info.information['channel']]['camera_type']
        self.config['sensor']=etc_info.information['cameras'][etc_info.information['channel']]['sensor']
        self.config['RN']=etc_info.information['cameras'][etc_info.information['channel']]['RN']
        self.config['DC']=etc_info.information['cameras'][etc_info.information['channel']]['DC']
        self.config['FWC']=etc_info.information['cameras'][etc_info.information['channel']]['FWC']
        self.config['gain']=etc_info.information['cameras'][etc_info.information['channel']]['gain']
        self.config['bits']=etc_info.information['cameras'][etc_info.information['channel']]['bits']
        self.config['xsize']=etc_info.information['cameras'][etc_info.information['channel']]['Nphotocell_X']
        self.config['ysize']=etc_info.information['cameras'][etc_info.information['channel']]['Nphotocell_Y']
        self.config['NrefPix_x']=etc_info.information['cameras'][etc_info.information['channel']]['NrefPix_x']
        self.config['NrefPix_y']=etc_info.information['cameras'][etc_info.information['channel']]['NrefPix_y']
        self.config['c1']=etc_info.information['cameras'][etc_info.information['channel']]['c1']
        self.config['c2']=etc_info.information['cameras'][etc_info.information['channel']]['c2']
        self.config['c3']=etc_info.information['cameras'][etc_info.information['channel']]['c3']
        self.config['c4']=etc_info.information['cameras'][etc_info.information['channel']]['c4']
        self.config['c5']=etc_info.information['cameras'][etc_info.information['channel']]['c5']
        self.config['c6']=etc_info.information['cameras'][etc_info.information['channel']]['c6']
        self.config['c7']=etc_info.information['cameras'][etc_info.information['channel']]['c7']
        self.config['c8']=etc_info.information['cameras'][etc_info.information['channel']]['c8']
        self.config['c9']=etc_info.information['cameras'][etc_info.information['channel']]['c9']
        self.config['xPixSize']=etc_info.information['cameras'][etc_info.information['channel']]['Photocell_SizeX']*etc_info.information['binning_X']
        self.config['yPixSize']=etc_info.information['cameras'][etc_info.information['channel']]['Photocell_SizeY']*etc_info.information['binning_Y']
        self.config['dig_noise']=etc_info.information['dig_noise']
        self.config['D_M1']=etc_info.information['D_M1']
        self.config['D_M2']=etc_info.information['D_M2']
        self.config['M2_factor']=etc_info.information['M2_factor']
        self.config['FoV_1axis']=etc_info.information['FoV_axis1']
        self.config['Focal_length']=etc_info.information['foc_len']
        self.config['filter_folder']=etc_info.information['filter_folder']
        self.config['Sky_CountRate']=etc_info.information['Sky_CountRate']
        #self.config['sky_brightness']=etc_info.information['sky_brightness']
        self.config['SB_eff']=etc_info.information['SB_eff']
        self.config['zeropoint']=etc_info.information['zeropoint']
        self.config['eff_wvl']=etc_info.information['effWavelength']
        self.config['pixelScale_X']=etc_info.information['pixelScale_X']
        self.config['pixelScale_Y']=etc_info.information['pixelScale_Y']
        self.config['airmass']=etc_info.information['airmass']
        self.config['seeing']=etc_info.information['seeing_los_arcsec']
        self.config['camera']=etc_info.information['channel']
        self.config['sky_site']=etc_info.information['sky_site']
        self.config['verbose']=str(etc_info.information['verbose'])

        if self.config['object_type'] == 'grb_sim':
            self.config['grb_mag']=etc_info.information['mag']

    def processConfigs(self):
        """
        Processes configuration information and save the information to a dictionary self.information.

        """
        #update the information dictionary
        self.information.update(self.config)

        #force gain to be float
        self.information['gain'] = float(self.config['gain'])

        # If resized image is present, need to change image size
        if 'ImageResized' in self.information:
            self.information['xsize']=self.information['ImageResized'][0]
            self.information['ysize']=self.information['ImageResized'][1]

        #name of the output file, include CCDs
        #self.information['output']=self.configfile[68:len(self.configfile)-6] + ".fits"

        #booleans to control the flow
        if self.config['shotNoise'].lower() == 'yes': self.shotNoise = True
        else: self.shotNoise = False
        if self.config['addSources'].lower() == 'yes': self.addsources = True
        else: self.addsources = False
        if self.config['bleeding'].lower() == 'yes': self.bleeding = True
        else: self.bleeding = False
        if self.config['cosmicRays'].lower() == 'yes': self.cosmicRays = True
        else: self.cosmicRays = False
        if self.config['cosmetics'].lower() == 'yes': self.cosmetics = True
        else: self.cosmetics = False

        #these don't need to be in the config file
        try:
            val=self.config['readoutNoise']
            if val.lower() == 'yes': self.readoutNoise = True
            else: self.readoutNoise = False
        except:
            self.readoutNoise = False
        try:
            val=self.config['digiNoise']
            if val.lower() == 'yes': self.digiNoise = True
            else: self.digiNoise = False
        except:
            self.digiNoise = False
        try:
            val=self.config['darkCurrent']
            if val.lower() == 'yes': self.darkCurrent = True
            else: self.darkCurrent = False
        except:
            self.darkCurrent = False
        try:
            val=self.config['persistance']
            if val.lower() == 'yes': self.persistance = True
            else: self.persistance = False
        except:
            self.persistance = False
        try:
            val=self.config['instrumIntrinsicNoise']
            if val.lower() == 'yes': self.instrumIntrinsicNoise = True
            else: self.instrumIntrinsicNoise = False
        except:
            self.instrumIntrinsicNoise = False
        try:
            val=self.config['background']
            if val.lower() == 'yes': self.background = True
            else: self.background = False
        except:
            self.background = False
        try:
            val=self.config['shutterOpen']
            if val.lower() == 'yes': self.shutterOpen = True
            else: self.shutterOpen = False
        except:
            self.shutterOpen = False
        try:
            val=self.config['nonLinearity']
            if val.lower() == 'yes': self.nonLinearity = True
            else: self.nonLinearity = False
        except:
            self.nonLinearity = False
        try:
            val=self.config['Vignetting']
            if val.lower() == 'yes': self.Vignetting = True
            else: self.Vignetting = False
        except:
            self.Vignetting = False
        try:
            val=self.config['ADU']
            if val.lower() == 'yes': self.ADU = True
            else: self.ADU = False
        except:
            self.ADU = False
        try:
            val=self.config['interpixCrosstalk']
            if val.lower() == 'yes': self.interpixCrosstalk = True
            else: self.interpixCrosstalk = False
        except:
            self.interpixCrosstalk = True
        try:
            val=self.config['Offset']
            if val.lower() == 'yes': self.Offset = True
            else: self.Offset = False
        except:
            self.Offset = False
        try:
            val=self.config['intscale']
            if val.lower() == 'yes': self.intscale = True
            else: self.intscale = False
        except:
            self.intscale = True

        self.information['variablePSF'] = False


        self.booleans = dict(shotNoiseoise=self.shotNoise,
                             addsources=self.addsources,
                             bleeding=self.bleeding,
                             cosmicRays=self.cosmicRays,
                             cosmetics=self.cosmetics,
                             background=self.background,
                             darkCurrent=self.darkCurrent,
                             readoutNoise=self.readoutNoise,
                             digiNoise=self.digiNoise,
                             nonLinearity=self.nonLinearity,
                             Vignetting=self.Vignetting,
                             ADU=self.ADU,
                             interpixCrosstalk=self.interpixCrosstalk,
                             persistance=self.persistance,
                             instrumIntrinsicNoise=self.instrumIntrinsicNoise,
                             Offset=self.Offset,
                             intscale=self.intscale,
                             shutterOpen=self.shutterOpen)

        if self.debug:
            pprint.pprint(self.information)

        
    def set_fits_header(self):
       """ Save information to save in FITS file header """
       from astropy.time import Time
       import datetime
       
       self.fits_header = OrderedDict()
       self.fits_header['FASTMODE'] = 0
       self.fits_header['NEXTRAP'] = 0
       self.fits_header['NEXTRAL'] = 0
       self.fits_header['REFOUT'] = 0
       self.fits_header['REFPIXEL'] = 0
       self.fits_header['ACQTIME'] = Time(datetime.datetime.now(), format='datetime', scale='utc').jd
       self.fits_header['ACQTIME1'] = str(datetime.datetime.now())
       self.fits_header['xsize'] = self.information['xsize']
       self.fits_header['xsize'] = self.information['xsize']
       self.fits_header['ysize'] = self.information['ysize']
       self.fits_header['FWC'] = self.information['FWC']
       self.fits_header['DC'] = round(self.information['DC'],3)
       self.fits_header['RN'] = self.information['RN']
       self.fits_header['gain'] = self.information['gain']
       self.fits_header['ZP'] = round(self.information['zeropoint'],3)
       self.fits_header['NEXP'] = self.information['Nexp']
       self.fits_header['EXPTIME'] = self.information['exptime']
       self.fits_header['SB'] = self.information['SB_eff']
       self.fits_header['airmass'] = round(self.information['airmass'],3)
       self.fits_header['seeing'] = round(self.information['seeing'],3)
       self.fits_header['SKYSITE'] = self.information['sky_site']
       self.fits_header['camera'] = self.information['camera']
       #self.fits_header['Temp_cam'] = self.information['Temp_cam']
       self.fits_header['bits'] = self.information['bits']
       self.fits_header['camera_type'] = self.information['camera_type']
       self.fits_header['sensor'] = self.information['sensor']
       self.fits_header['binning_X'] = self.information['binning_X']
       self.fits_header['binning_Y'] = self.information['binning_Y']
       self.fits_header['filter'] = self.information['filter_folder']+'-'+self.information['filter_band']
       self.fits_header['D_M1'] = self.information['D_M1']
       self.fits_header['D_M2'] = self.information['D_M2']
       self.fits_header['M2_factor'] = self.information['M2_factor']
       if self.acquisition=='ramp' : 
           self.fits_header['GROUP'] = self.nbGroup
           self.fits_header['READ'] = self.nbRead
           self.fits_header['DROP'] = self.nbDrop

       # add WCS to the header
       self.fits_header['WCSAXES'] = 2
       self.fits_header['CRPIX1'] = self.information['ysize']/2.
       self.fits_header['CRPIX2'] = self.information['xsize']/2.
       self.fits_header['CRVAL1'] = self.information['RA']
       self.fits_header['CRVAL2'] = self.information['DEC']
       self.fits_header['CTYPE1'] = 'RA---TAN'
       self.fits_header['CTYPE2'] = 'DEC--TAN'
       #north is up, east is left
       #self.fits_header['CD1_1'] = -self.config['pixelScale_X'] / 3600. #pix size in arc sec / deg
       #self.fits_header['CD1_2'] = 0.0
       #self.fits_header['CD2_1'] = 0.0
       #self.fits_header['CD2_2'] = self.config['pixelScale_Y'] / 3600.

       self.fits_header['CDELT1'] = -self.config['pixelScale_X'] / 3600. #pix size in arc sec / deg
       self.fits_header['CDELT2'] = self.config['pixelScale_Y'] / 3600.
       self.fits_header['CROTA2'] = 0.0
       self.fits_header['DATE-OBS'] = datetime.datetime.isoformat(datetime.datetime.now())
       self.fits_header['INSTRUME'] = 'ImSimpy' 


       #create a new FITS file, using HDUList instance
       hdu=fits.PrimaryHDU()

       #new image HDU
       #hdu.data=self.image.astype(np.float32)
 
       #add input keywords to the header
       for key, value in self.fits_header.items():
           #truncate long keys
           if len(key) > 8:
               key = key[:7]
           try:
               hdu.header.set(key.upper(), value)
           except:
               try:
                   hdu.header.set(key.upper(), str(value))
               except:
                   pass

       #write booleans
       for key, value in self.booleans.items():
           #truncate long keys
           if len(key) > 8:
               key = key[:7]
           hdu.header.set(key.upper(), str(value), 'Boolean Flags')

       hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
       hdu.header.add_history('Created by ImSimpy at %s' % datetime.datetime.isoformat(datetime.datetime.now()))
       hdu.verify('fix')

       # Create directory if not existing
       os.makedirs(os.path.dirname(self.path+'/images/'+self.information['output']),exist_ok=True)
       
       #write the actual file
       #hdu.writeto(self.path +'/images/'+self.information['output'],overwrite=True)
       #print (hdu.header)
       self.hdu_header=hdu.header


    def _createEmpty(self):
        """
        Creates and empty array of a given x and y size full of zeros.
        """
        self.image_total = np.zeros((self.information['ysize'], self.information['xsize']), dtype=np.float64)
       
        #reference pixel case management
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            self.image = self.image_total[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
        else:
            self.image = self.image_total


    def objectOnDetector(self, object):
        """
        Tests if the object falls on the detector area being simulated.

        :param object: object to be placed to the self.image being simulated.
        :type object: list

        :return: whether the object falls on the detector or not
        :rtype: bool
        """
        ny, nx = self.finemap[object[3]].shape
        xt = object[0]
        yt = object[1]

        #the bounding box should be in the nominal scale
        fac = 1./self.information['psfoversampling']

        #Assess the boundary box of the input image
        xlo = (1 - nx) * 0.5 * fac + xt
        xhi = (nx - 1) * 0.5 * fac + xt
        ylo = (1 - ny) * 0.5 * fac + yt
        yhi = (ny - 1) * 0.5 * fac + yt

        i1 = np.floor(xlo + 0.5)
        i2 = np.ceil(xhi + 0.5) + 1
        j1 = np.floor(ylo + 0.5)
        j2 = np.ceil(yhi + 0.5) + 1

        if i2 < 1 or i1 > self.information['xsize']:
            return False

        if j2 < 1 or j1 > self.information['ysize']:
            return False

        return True


    def overlayToCCD(self, data, obj):
        """
        Overlay data from a source object onto the self.image

        :param data: ndarray of data to be overlaid on to self.image
        :type data: ndarray
        :param obj: object information such as x,y position
        :type obj: list
        """
        #object centre x and y coordinates (only in full pixels, fractional has been taken into account already)
        xt = np.floor(obj[0]) - 1  #zero indexing
        yt = np.floor(obj[1]) - 1  #zero indexing

        #input array size
        nx = data.shape[1]
        ny = data.shape[0]

        # Assess the boundary box of the input image
        xlo = (1 - nx) * 0.5 + xt
        xhi = (nx - 1) * 0.5 + xt + 1
        ylo = (1 - ny) * 0.5 + yt
        yhi = (ny - 1) * 0.5 + yt + 1

        i1 = int(np.floor(xlo + 0.5))
        if i1 < 1:
            i1 = 0

        i2 = int(np.floor(xhi + 0.5))
        if i2 > self.information['xsize']:
            i2 = self.information['xsize']

        j1 = int(np.floor(ylo + 0.5))
        if j1 < 1:
            j1 = 0

        j2 = int(np.floor(yhi + 0.5))
        if j2 > self.information['ysize']:
            j2 = self.information['ysize']

        if i1 > i2 or j1 > j2:
            #print ('Object does not fall on the detector...')
            return

        ni = i2 - i1
        nj = j2 - j1

        
        #add to the image
        if ni == nx and nj == ny:
            #full frame will fit
            self.image[j1:j2, i1:i2] += data
        elif ni < nx and nj == ny:
            #x dimensions shorter
            if int(np.floor(xlo + 0.5)) < 1:
                #small values, left side
                self.image[j1:j2, i1:i2] += data[:, nx-ni:]
            else:
                #large values, right side
                self.image[j1:j2, i1:i2] += data[:, :ni]
        elif nj < ny and ni == nx:
            #y dimensions shorter
            if int(np.floor(ylo + 0.5)) < 1:
                #small values, bottom
                self.image[j1:j2, i1:i2] += data[ny-nj:, :]
            else:
                #large values, top
                self.image[j1:j2, i1:i2] += data[:nj, :]
        else:
            #both lengths smaller, can be in any of the four corners
            if int(np.floor(xlo + 0.5)) < 1 > int(np.floor(ylo + 0.5)):
                #left lower
                self.image[j1:j2, i1:i2] += data[ny-nj:, nx-ni:]
            elif int(np.floor(xlo + 0.5)) < 1 and int(np.floor(yhi + 0.5)) > self.information['ysize']:
                #left upper
                self.image[j1:j2, i1:i2] += data[:nj, nx-ni:]
            elif int(np.floor(xhi + 0.5)) > self.information['xsize'] and int(np.floor(ylo + 0.5)) < 1:
                #right lower
                self.image[j1:j2, i1:i2] += data[ny-nj:, :ni]
            else:
                #right upper
                self.image[j1:j2, i1:i2] += data[:nj, :ni]


    def writeFITSfile(self, data, filename, unsigned16bit=False):
        """
        Writes out a simple FITS file.

        :param data: data to be written
        :type data: ndarray
        :param filename: name of the output file
        :type filename: str
        :param unsigned16bit: whether to scale the data using bzero=32768
        :type unsigned16bit: bool

        :return: None
        """
        if os.path.isfile(filename):
            os.remove(filename)

        #create a new FITS file, using HDUList instance
        #hdulist = fits.HDUList(fits.PrimaryHDU())
        hdu=fits.PrimaryHDU()
 

        #new image HDU
        #hdu = fits.ImageHDU(data=data)
        hdu.data=data

        #convert to unsigned 16bit int if requested
        if unsigned16bit:
            hdu.scale('int16', '', bzero=32768)
            hdu.header.add_history('Scaled to unsigned 16bit integer!')

        #add input keywords to the header
        for key, value in self.fits_header.items():
            #truncate long keys
            if len(key) > 8:
                key = key[:7]
            try:
                hdu.header.set(key.upper(), value)
            except:
                try:
                    hdu.header.set(key.upper(), str(value))
                except:
                    pass

        #write booleans
        for key, value in self.booleans.items():
            #truncate long keys
            if len(key) > 8:
                key = key[:7]
            hdu.header.set(key.upper(), str(value), 'Boolean Flags')

        #update and verify the header
        hdu.header.add_history('This is an intermediate data product no the final output!')
        hdu.header.add_history('Created by ImSimpy at %s' % (                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
        hdu.verify('fix')

        #hdulist.append(hdu)

        #write the actual file
        #hdulist.writeto(filename)
        hdu.writeto(filename,overwrite=True)


    def configure(self,config_type):
        """
        Configures the simulator with input information and creates and empty array to which the final image will
        be build on.
        """
        if config_type=='file': self.readConfigs()
        self.etc(config_type)
        self.processConfigs()
        self._createEmpty()
        self.set_fits_header()


    def generateObjectList(self):
       """ Generate object to simulate """

       if 'generate' in self.config['SourcesList']:
           from astropy.io import fits

           if 'output' in self.information['SourcesList']['generate']:
               output=self.path+'/data/catalog/'+self.information['SourcesList']['generate']['output']
           else:
               output=self.path+'/data/catalog/SourcesCatalog.txt'
           if 'frame' in self.information['SourcesList']['generate']:
               frame=self.information['SourcesList']['generate']['frame']
           else:
               frame='icrs'
           if 'band' in self.information['SourcesList']['generate']:
               band = self.information['SourcesList']['generate']['band']
           else:
               band=self.information['filter_band']

           RA=self.information['SourcesList']['generate']['RA']
           DEC=self.information['SourcesList']['generate']['DEC']
           radius=self.information['SourcesList']['generate']['radius']
           #_header=fits.open(self.path+'/images/'+self.information['output'])
           #header=_header['PRIMARY'].header

           if self.information['SourcesList']['generate']['catalog'] == 'Panstarrs':
               print ('Downloading objects from Panstarrs catalog')
               from ImSimpy.utils.createCatalogue import  PanstarrsCatalog
               PanstarrsCatalog(RA, DEC, radius, band, self.config['eff_wvl'], self.hdu_header, frame=frame, output=output)
           else:
               from ImSimpy.utils.createCatalogue import  Viziercatalog
               print ('Downloading objects from Vizier')
               Viziercatalog(RA, DEC, radius, band, self.config['eff_wvl'], self.hdu_header, catalog=self.information['SourcesList']['generate']['catalog'], frame=frame, output=output)
           self.objects = np.loadtxt(output)

       elif "file" in self.config['SourcesList']:
           self.objects = np.loadtxt(self.path+'/data/catalog/'+self.information['SourcesList']['file'])


    def readObjectlist(self):
        """
        Reads object list using numpy.loadtxt, determines the number of object types,
        and finds the file that corresponds to a given object type.

        The input catalog is assumed to contain the following columns:

            #. x coordinate
            #. y coordinate
            #. apparent magnitude of the object
            #. type of the object [0=star, number=type defined in the objects.dat]
            #. rotation [0 for stars, [0, 360] for galaxies]

        This method also displaces the object coordinates based on the
        CCD to be simulated.

        .. Note:: If even a single object type does not have a corresponding input then this method
                  forces the program to exit.
        """
        #self.objects = np.loadtxt(self.path+self.information['SourcesList'])
        #Add GRB on the object list
        if self.config['object_type'] == 'grb_sim':
           # If coordinates given in pixels
           if self.config['grb_coord_type'] == 'pixels':
               grb_coord_pix=self.config['grb_coords']
           elif self.config['grb_coord_type'] == 'RADEC':
               from astropy.io import fits
               import astropy.units as u
               import astropy.coordinates as coord
               from astropy.wcs import WCS
               
               c = coord.SkyCoord(self.config['grb_coords'][0], self.config['grb_coords'][1], unit=(u.deg, u.deg),frame='icrs')
               #_header=fits.open(self.path + '/images/' +self.information['output'])
               #header=_header['PRIMARY'].header
               
               w = WCS(self.hdu_header)
               world = np.array([[c.ra.deg, c.dec.deg]])
               #print (world)
               pix = w.all_world2pix(world,1)
               #print (pix) 
               #first transform WCS into pixels
               grb_coord_pix=pix[0]
           self.config['grb_coords_pix_X']=grb_coord_pix[1]
           self.config['grb_coords_pix_Y']=grb_coord_pix[0]
           self.objects = np.vstack((self.objects, [grb_coord_pix[0],grb_coord_pix[1],self.config['grb_mag'],1000,0]))
           #Add GRB to the object list as a point source
           txt = 'GRB positionned at pixel coordinates (X,Y): {0:.2f},{1:.2f} with mag= {2:.2f}'.format(grb_coord_pix[0],grb_coord_pix[1],round(self.config['grb_mag'],2))
           print (txt)

              
        #if only a single object in the input, must force it to 2D
        try:
            tmp_ = self.objects.shape[1]
        except:
            self.objects = self.objects[np.newaxis, :]

        #read in object types
        data = open(self.path+'/data/objects.dat').readlines()

        #only 2D array will have second dimension, so this will trigger the exception if only one input source
        tmp_ = self.objects.shape[1]
        #find all object types
        self.sp = np.asarray(np.unique(self.objects[:, 3]), dtype=np.int)

        #generate mapping between object type and data
        objectMapping = {}
        for stype in self.sp:
            
            if stype == 0 or stype == 1000:
                #delta function
                objectMapping[stype] = 'PSF'
            else:
                for line in data:
                    tmp = line.split()
                    if int(tmp[0]) == stype:
                        #found match
                        if tmp[2].endswith('.fits'):
                            d,header = fits.getdata(self.path+'/'+tmp[2],header=True)
                            #print (type(d),d.shape,d,np.max(d))
                            if 'PIXSIZE1' in header:
                                ratio=float(header['PIXSIZE1']/(self.information['xPixSize']*1e6))
                                print ('ratio finemaps: %.2f' % ratio)
                                d2=scipy.ndimage.zoom(d/d.sum(), ratio, order=3)
                                d2[d2<1e-6]=0
                                if np.sum(d2) != 1: d2=d2/np.sum(d2)
                                image_size=d2.shape
                                #Assume same size horizontally and vertically
                                if image_size[0] % 2 == 0:
                                    width=int(image_size[0]/ratio/4)
                                else:
                                    width=int((image_size[0]/ratio-1)/4)

                                #Assume same size horizontally and vertically
                                if np.ceil(image_size[0]) % 2 == 0:
                                    center=int((np.ceil(image_size[0]))/2)
                                else:
                                    center=int((np.ceil(image_size[0])-1)/2)
                                center=[center,center]
                                #print (ratio)
                                #print (center)
                                #print (width)
                                #print (center[0]-width,center[0]+width,center[1]-width,center[1]+width)
                                d3=d2[center[0]-width:center[0]+width,center[1]-width:center[1]+width]
                                

                                #print (type(d3),d3.shape,d3,np.max(d3))
                            else:
                                print ('No pixel size found in header. Assume the same as the current telescope.')

                                d3=d
                        else:
                            pass
                        objectMapping[stype] = dict(file=tmp[2], data=d3)
                        break

        self.objectMapping = objectMapping


    def generatePSF(self):
       """ Compute PSF if needed """

       PSF = dict()
       # Load atmosphere and instrument PSF 
       if self.information['PSF']['total']['method'] == 'compute': 
           for keyword in self.information['PSF']:
               if keyword != "total":
                   if "file" not in self.information['PSF'][keyword]:
                       if self.information['PSF'][keyword]['type'] == 'moffat':
                           if 'beta' in self.information['PSF'][keyword]:
                               beta=self.information['PSF'][keyword]['beta']
                           else:
                               beta=2
                       else:
                           beta=2
                       if 'seeing' in self.information['PSF'][keyword]:
                               seeing=self.information['PSF'][keyword]['seeing']
                       else:
                               seeing=self.config['seeing']

                       # If PSF size bigger than image --> Limit PSF size to image size
                       if self.information['PSF'][keyword]['size'][0] > self.information['xsize']:
                           self.information['PSF'][keyword]['size'][0] = self.information['xsize']
                           print ('PSF size along x axis bigger than image size!\nPSF size limited to image size along x axis now: %d Pixels' % (self.information['xsize']))

                       if self.information['PSF'][keyword]['size'][1] > self.information['ysize']:
                           self.information['PSF'][keyword]['size'][1] = self.information['ysize']
                           print ('PSF size along y axis bigger than image size!\nPSF size limited to image size along y axis now: %d Pixels' % (self.information['ysize']))

                       PSFUtils.createPSF(filename=self.path+'/data/psf/'+self.information['PSF'][keyword]['output'],PSF_type=self.information['PSF'][keyword]['type'],imsize=self.information['PSF'][keyword]['size'],pixel_size=[self.config['xPixSize'],self.config['yPixSize']],pixel_scale=self.config['pixelScale_X'],eff_wvl=self.config['eff_wvl'],seeing=seeing,DM1=self.config['D_M1'],DM2=self.config['D_M2'],focal_length=self.config['Focal_length'],oversamp=self.config['psfoversampling'],beta=beta,disp=False,unsigned16bit=False)
             
                       PSF[keyword] = self.path+'/data/psf/'+self.information['PSF'][keyword]['output']

                   else:
                       # Check pixel size and oversample if needed
                       hdr_ = fits.getheader(self.path+'/data/psf/'+self.information['PSF'][keyword]['file']+'.fits')
                       try: 
                           if hdr_['XPIXELSZ'] != self.information['cameras'][self.information['channel']]['Photocell_SizeX'] / oversamp or hdr_['YPIXELSZ'] != self.information['cameras'][self.information['channel']]['Photocell_SizeY'] / oversamp :
                               resampling=[self.information['cameras'][self.information['channel']]['Photocell_SizeX'] / oversamp,self.information['cameras'][self.information['channel']]['Photocell_SizeY'] / oversamp]

                               PSFUtils.resize(filename1=self.path+'/data/psf/'+self.information['PSF']['keyword']['file'],filename2=self.path+self.information['PSF']['keyword']['file']+'_oversammpled',type='factor',resampling=resampling,overwrite=True,unsigned16bit=False)

                               PSF[keyword] = self.path+'/data/psf/'+self.information['PSF'][keyword]['file']+'_oversampled'
                           else:
                               PSF[keyword] = self.path+'/data/psf/'+self.information['PSF'][keyword]['file']
                       except:
                               PSF[keyword] = self.path+'/data/psf/'+self.information['PSF'][keyword]['file']
           print ('PSF convolution')
           # convolve atmosphere and instrument PSF to get the total PSF
           PSFUtils.convolvePSF(filename1=PSF['atmosphere'],filename2=PSF['instrument'],filename3=self.path+'/data/psf/'+self.information['PSF']['total']['file'])
           #PSFUtils.convolvePSF(filename1=PSF['instrument'],filename2=PSF['atmosphere'],filename3=self.path+self.information['PSF']['total']['output']+'_oversampled')
           #PSFUtils.resize(filename1=self.path+self.information['PSF']['total']['output']+'_oversampled',filename2=self.path+self.information['PSF']['total']['output'],resampling=32/self.information['psfoversampling'],type='sum')
           #PSFUtils.resize(filename1=self.path+self.information['PSF']['total']['output']+'_oversampled',filename2=self.path+self.information['PSF']['total']['output'],resampling=self.information['psfoversampling']/32,type='zoom')
           print ('done')

    def readPSFs(self):
        """
        Reads in a PSF from a FITS file.

        .. Note:: at the moment this method supports only a single PSF file.
        """
        if self.information['variablePSF']:
            #grid of PSFs
            print ('Spatially variable PSF not implemented -- exiting')
            sys.exit(-9)
        else:
            #single PSF
            self.PSF = fits.getdata(self.path+'/data/psf/'+self.information['PSF']['total']['file']).astype(np.float64)
            # Normalise if needed
            if np.sum(self.PSF) != 1 : self.PSF /= np.sum(self.PSF)
            self.PSFx = self.PSF.shape[1]
            self.PSFy = self.PSF.shape[0]


    def generateFinemaps(self):
        """
        Generates finely sampled images of the input data.
        """
        self.finemap = {}
        self.shapex = {}
        self.shapey = {}
        for k, stype in enumerate(self.sp):
            if stype == 0 or stype == 1000:
                data = self.PSF.copy().astype(np.float64)
                data /= np.sum(data)
                self.finemap[stype] = data
                self.shapex[stype] = 0
                self.shapey[stype] = 0
            else:

                # Rescaled to pixel size
                


                if self.information['psfoversampling'] > 1.0:
                    data = scipy.ndimage.zoom(self.objectMapping[stype]['data'],
                                              self.information['psfoversampling'],
                                              order=0)
                else:
                    data = self.objectMapping[stype]['data']

                data[data < 0.] = 0.0
                if data.sum() != 1: data /= np.sum(data)
                self.finemap[stype] = data



    def addObjects(self):
        """
        Add objects from the object list to the CCD image (self.image).

        Scale the object's brightness in electrons and size using the input catalog magnitude.
        
        """
        #total number of objects in the input catalogue and counter for visible objects
        n_objects = self.objects.shape[0]
        visible = 0

        print ('Total number of objects in the input catalog = %i' % n_objects)

        #calculate the scaling factors from the magnitudes
        #intscales = 10.0**(-0.4 * self.objects[:, 2])*self.information['magzero']) * self.information['exptime']
        #calculate the number of electrons from the magnitudes
        mag2elec = 10.0**(-0.4 * (self.objects[:, 2] -  self.information['zeropoint'])) * self.information['exptime']

        intscales=mag2elec
        if ~self.random:
            #Using a fixed size-magnitude relation (equation B1 from Miller et al. 2012 (1210.8201v1).

            #testin mode will bypass the small random scaling in the size-mag relation
            #loop over exposures
            for i in range(self.information['Nexp']):
                #loop over the number of objects
                for j, obj in enumerate(self.objects):

                    stype = obj[3]

                    if self.objectOnDetector(obj):
                        visible += 1
                        if stype == 0 or stype == 1000:
                            #point source, apply PSF
                            if stype == 0:
                                txt = "Star: " + str(j + 1) + "/" + str(n_objects) + \
                                  " mag=" +str(round(obj[2],2)) + " nb of el.=" + str(round(mag2elec[j],2))
                            elif stype == 1000:
                                txt = "GRB: " + str(j + 1) + "/" + str(n_objects) + \
                                  " mag=" +str(round(obj[2],2)) + " nb of el.=" + str(round(mag2elec[j],2))

                            if self.debug: print (txt)

                            data = self.finemap[stype].copy()
                            #print (data.shape) 
                            #map the data to new grid aligned with the centre of the object and scale
                            yind, xind = np.indices(data.shape)
                            if self.information['psfoversampling'] != 1.0:
                                yi = yind.astype(np.float) + self.information['psfoversampling']*(obj[0] % 1)
                                xi = xind.astype(np.float) + self.information['psfoversampling']*(obj[1] % 1)
                            else:
                                yi = yind.astype(np.float) + (obj[0] % 1)
                                xi = xind.astype(np.float) + (obj[1] % 1)
                            data = ndimage.map_coordinates(data, [yi, xi], order=1, mode='nearest')
                            if self.information['psfoversampling'] != 1.0:
                                data = scipy.ndimage.zoom(data, 1. / self.information['psfoversampling'], order=1)
                            #suppress negative numbers, renormalise and scale with the intscale
                            data[data < 0.0] = 0.0
                            sum = np.sum(data)
                            sca = mag2elec[j] / sum

                            """
                            print ('Obj coord: %.2f %.2f  / mag: %.2f   / finemap: min: %f  max: %f  mean: %f' % (obj[0],obj[1],obj[2],np.min(data*sca),np.max(data*sca),np.mean(data*sca)))
                            tolerance_min_finemap=1e-2
                            test_min=np.min(data*sca)
                            zoom=2
                            while test_min > tolerance_min_finemap:
                               print ('increased fineap map by 2 ')
                               data = scipy.ndimage.zoom(data, zoom)
                               test_min=np.min(data*sca)
                               print (np.min(data*sca))
                            """
                            #sca = intscales[j] / sum
                            #data = ne.evaluate("data * sca")
                            #sca = mag2elec[j]

                            # numexpr apparently faster than numpy for big arrays
                            data = ne.evaluate("data * sca")

                            #print ('Obj coord: %.2f %.2f  / mag: %.2f   / finemap: min: %f  max: %f  mean: %f' % (obj[0],obj[1],obj[2],np.min(data),np.max(data),np.mean(data)))
                            #data[data < 0.0] = 0.0

                            #overlay the scaled PSF on the image
                            self.overlayToCCD(data, obj)
                        else:
                            #extended source, rename finemap
                            data = self.finemap[stype].copy()
                            #map the data to new grid aligned with the centre of the object
                            yind, xind = np.indices(data.shape)
                            if self.information['psfoversampling'] != 1.0:
                                yi = yind.astype(np.float) + self.information['psfoversampling']*(obj[0] % 1)
                                xi = xind.astype(np.float) + self.information['psfoversampling']*(obj[1] % 1)
                            else:
                                yi = yind.astype(np.float) + (obj[0] % 1)
                                xi = xind.astype(np.float) + (obj[1] % 1)

                            #yi = yind.astype(np.float) + (obj[0] % 1)
                            #xi = xind.astype(np.float) + (obj[1] % 1)
                            data = ndimage.map_coordinates(data, [yi, xi], order=1, mode='nearest')

                            conv = convolution(data, self.PSF, self.information['mode'])

                            #suppress negative numbers
                            conv[conv < 0.0] = 0.0

                            #renormalise and scale to the right magnitude
                            sum = np.sum(conv)
                            #sca = intscales[j] / sum
                            sca = mag2elec[j]/sum
                            conv = ne.evaluate("conv * sca")

                            #tiny galaxies sometimes end up with completely zero array
                            #checking this costs time, so perhaps this could be removed
                            if np.isnan(np.sum(conv)):
                                continue

                            #overlay the convolved image on the image
                            self.overlayToCCD(conv, obj)
                            
                    else:
                        #not on the screen
                        #print ('Object %i is outside the detector area' % (j + 1))
                        pass
                        
            print ('%i/%i objects were placed on the detector' % (visible,n_objects))



    def addCosmetics(self):
        """ Add the cosmetics effects """
        deadPixs=fits.getdata(self.path+'/data/Cosmetics/'+self.information['DeadPixFile'])
        HotPixs=fits.getdata(self.path+'/data/Cosmetics/'+self.information['HotPixFile'])
        self.image*=deadPixs
        self.image+=HotPixs


    def addCosmicRays(self):
        """ Add cosmic rays """
        pass


    def applyVignetting(self):
        """ Add vignetting  """
        vignetting=fits.getdata(self.path+'/data/Vignetting/Calibration/colibri')
        self.image*=vignetting


    def applyNonLinearity(self):
        """ Add non linearity  """
        import copy
        A1=fits.getdata(self.path+'\\data\\NonLinearity\\'+self.information['nonLinearity_A1'])
        A2=fits.getdata(self.path+'\\data\\NonLinearity\\'+self.information['nonLinearity_A2'])
        gain_conv=fits.getdata(self.path+'/data/GainMap/'+self.information['GainMapFile'])

        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            A1 = A1[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
            A2 = A2[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
            gain_conv = gain_conv[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
        A2[A2 > 0.0]=0.0
        A1[A1==0.0]=0.01
        
        image=copy.deepcopy(self.image)
        self.image+=A2/((A1**2)*gain_conv)*image**2
        
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            self.image_total[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]=self.image

        

    def applyDarkCurrent(self):
        """
        Apply dark current. Scales the dark with the exposure time. Apply on reference pixel.
        """
        filename_DC=self.path+'\\data\\DarkCurrent\\'+'DarkCurrent_'+self.information['camera']+'.fits'
        DarkCurrent(filename=filename_DC,Type='gaussian',mean=self.information['DC'],xsize=self.information['xsize'],ysize=self.information['ysize'],texp=self.information['exptime'])

        DC=fits.getdata(filename_DC)

        self.image_total+=DC


    def applySkyBackground(self):
        """
        Apply dark the sky background. Scales the background with the exposure time.

        Additionally saves the image without noise to a FITS file.
        """
        #np.random.seed(seed=1234567890)

        #add background
        #sky_back_el = 10**(-0.4*(self.information['sky_brightness']-self.information['zeropoint']))*self.information['pixelScale']**2
        sky_back_el=self.information['Sky_CountRate']
        #print (sky_back_el, self.information['Sky_CountRate'])
        #bcgr = self.information['exptime'] * self.information['cosmic_bkgd']
        bcgr = self.information['exptime'] * sky_back_el
        #self.image += bcgr + np.sqrt(np.random.poisson(bcgr),size=self.image.shape)
        sky_image = np.random.poisson(bcgr,size=self.image.shape).astype(np.float64)
        self.image += sky_image
        
        
    def applyRampSkyBackground(self):
        """
        Apply dark the sky background for a ramp. Scales the background with the exposure time.

        Additionally saves the image without noise to a FITS file.
        """
        sky_back_el=self.information['Sky_CountRate']
        bcgr = 1.475 * sky_back_el
        sky_image = np.random.poisson(bcgr,size=self.image.shape).astype(np.float64)
        self.image += sky_image

        

    def applyShotNoise(self):
        """
        Add Poisson noise to the image.
        """

        self.image = np.random.poisson(self.image).astype(np.float64)




    def applyReadoutNoise(self):
        """
        Applies readout noise to the image being constructed.

        The noise is drawn from a Normal (Gaussian) distribution with average=0.0 and std=readout noise.
        """

        noise = np.random.normal(loc=0.0, scale=self.information['RN'], size=self.image.shape)
        #noise = np.random.poisson(loc=0.0,scale=self.information['RN'], size=self.image.shape)

        #add to the image
        self.image += noise
        
        # Replace self.image in self.image_total
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            self.image_total[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]=self.image


    def applyDigiNoise(self):
        """
        Applies digitisation noise to the image being constructed.

        """
        #noise = np.random.normal(loc=0.0, scale=self.information['RN'], size=self.image.shape)
        diginoise = np.random.poisson(self.information['gain']/np.sqrt(12), size=self.image.shape)

        # Can not be negative
        diginoise[diginoise < 0.0] = 0.0

        #add to the image
        self.image += diginoise


    def electrons2ADU(self):
        """
        Convert from electrons to ADUs using the value read from the configuration file.
        """
        gain_map=fits.getdata(self.path+'/data/GainMap/'+self.information['GainMapFile'])
        
        self.image_total /= gain_map
        
        # Replace self.image in self.image_total
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            self.image=self.image_total[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]


    def addOffset(self):
        """
        Add the offset (bias) in ADU
        """
        offset=fits.getdata(self.path+'/data/Offset/'+self.information['OffsetFile'])
        gain_conv=fits.getdata(self.path+'/data/GainMap/'+self.information['GainMapFile'])
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            offset = offset[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
            gain_conv=gain_conv[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
        
        self.image+=offset*gain_conv

    
    def applyBleeding(self):
        """
        Apply bleeding along the CCD columns if the number of electrons in a pixel exceeds the full-well capacity.

        Bleeding is modelled in the parallel direction only, because the CCD273s are assumed not to bleed in
        serial direction.

        :return: None
        """
        if self.config['verbose'] == 'True': 
            print ('Max el. in 1 pix: %f, FWC: %i ' % (np.max(self.image),self.information['FWC']))
        if np.any(self.image > self.information['FWC']):
        
            #loop over each column, as bleeding is modelled column-wise
            for i, column in enumerate(self.image.T):
                sum = 0.
                for j, value in enumerate(column):
                    #first round - from bottom to top (need to half the bleeding)
                    overload = value - self.information['FWC']
                    if overload > 0.:
                        overload /= 2.
                        self.image[j, i] -= overload
                        sum += overload
                    elif sum > 0.:
                        if -overload > sum:
                            overload = -sum
                        self.image[j, i] -= overload
                        sum += overload

            for i, column in enumerate(self.image.T):
                sum = 0.
                for j, value in enumerate(column[::-1]):
                    #second round - from top to bottom (bleeding was half'd already, so now full)
                    overload = value - self.information['FWC']
                    if overload > 0.:
                        self.image[-j-1, i] -= overload
                        sum += overload
                    elif sum > 0.:
                        if -overload > sum:
                            overload = -sum
                        self.image[-j-1, i] -= overload
                        sum += overload


    def applyInterpixCrosstalk(self):
        """
        Apply interpixel crosstalk coefficients to the image with a convolution
        """
        import copy
        ICTcoeff = np.array([[self.information['c1'], self.information['c2'],self.information['c3']],[self.information['c4'],self.information['c5'],self.information['c6']],[self.information['c7'],self.information['c8'],self.information['c9']]])
        
        # apply coefficients to the image
        copie=np.float32(copy.deepcopy(self.image))
        convolve = ndimage.convolve(copie, ICTcoeff, mode='constant', cval=0.0)   
        self.image = copy.deepcopy(convolve)    
        
        # Replace self.image in self.image_total
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            self.image_total[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]=self.image
        else :
            self.image_total=self.image

        
    def applyInstrumIntrinsicNoise(self) :
        """
        Apply the instrument intrinsic noise
        """
        instrumIntrinsicNoise=fits.getdata(self.path+'/data/IntrinsicNoise/'+self.information['instrumIntrinsicNoiseFile'])
        gain_conv=fits.getdata(self.path+'/data/GainMap/'+self.information['GainMapFile'])
        #reference pixel case management
        if 'NrefPix_x' in self.config and 'NrefPix_y' in self.config :
            instrumIntrinsicNoise = instrumIntrinsicNoise[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
            gain_conv=gain_conv[self.information['NrefPix_x']:-self.information['NrefPix_x'], self.information['NrefPix_y']:-self.information['NrefPix_y']]
        
        self.image+=instrumIntrinsicNoise*self.config['exptime']*gain_conv

    
    def applyPersistance(self):
        """Apply persistance"""
        pass


    def discretise(self):
        """
        Converts a floating point image array (self.image) to an integer array with max values
        defined by the argument max.

        :param max: maximum value the the integer array may contain [default 65k]
        :type max: float

        :return: None
        """
        max=2**self.information['bits']-1

        #cut of the values larger than max
        self.image[self.image > max] = max

        self.image = np.rint(self.image).astype(np.int)


    def writeOutputs(self):
        """
        Writes out a FITS file using PyFITS and converts the image array to 16bit unsigned integer as
        appropriate for VIS.

        Updates header with the input values and flags used during simulation.
        """
        #if os.path.isfile(self.path +  self.information['output_dir']+self.information['output']):
        #    os.remove( self.path + self.information['output_dir']+self.information['output'])
 
        
        #create a new FITS file, using HDUList instance
        hdu=fits.PrimaryHDU()
        
        #hdu=fits.open(self.path + '/images/'+self.information['output'])

        #new image HDU
        #float 32 bits
        #hdu['PRIMARY'].data=self.image.astype(np.float32)
        # UNsigned Integer 16 bits (0 to 65535)
        hdu.data=self.image_total.astype(np.uint16)
        hdu.header=self.hdu_header

        #write the actual file. Path should exists, already checked while creating the headers.
        hdu.writeto(self.path + '/images/'+self.information['output'],overwrite=True)


    def simulate(self,config_type='file'):
        """
        Create a single simulated image defined by the configuration file.
        Will do all steps defined in the config file sequentially.

        :return: None
        """
        #if self.config['verbose'] == 'True': print ("Read config file and execute ETC")
        print ("Read config file and execute ETC")
        #self.information['output']=
        self.configure(config_type)
        self.information['output']=self.nom + ".fits"
        print ("Building image: %s:" % self.information['output'])
        #print (self.information)
               
        if self.Offset:
            #if self.config['verbose'] == 'True': print ("Add Offset")
            print ("\tAdd Offset")
            self.addOffset()
        
        if self.addsources:
            #if self.config['verbose'] == 'True': print ("Read objecct list")
            print ("\tGENERATE OBJECTS CATALOG")
            self.generateObjectList()
            self.readObjectlist()
        
            #if self.config['verbose'] == 'True': print ("Generate PSF")
            print ("\tGENERATE PSF")
            self.generatePSF()
            self.readPSFs()
            self.generateFinemaps()
       
            #if self.config['verbose'] == 'True': print ("Add objects")
            print ("\tADD OBJECTS")
            self.addObjects()

        if self.background:
            #if self.config['verbose'] == 'True': print ("Add Sky background")
            print ("\tAdd Sky background")
            self.applySkyBackground()

        if self.Vignetting:
            #if self.config['verbose'] == 'True': print ("Add Vignetting")
            print ("\tAdd Vignetting")
            self.applyVignetting()

        if self.cosmicRays:
            #if self.config['verbose'] == 'True': print ("Add cosmic Rays")
            print ("\tAdd cosmic Rays")
            self.addCosmicRays()
        """ N'apparait pas, voir tuto
        if self.shotNoise:
            #if self.config['verbose'] == 'True': print ("Apply Shot noise")
            print ("\tApply Shot noise")
            self.applyShotNoise()
        """
        if self.darkCurrent:
            #if self.config['verbose'] == 'True': print ("Add dark current")
            print ("\tAdd dark current")
            self.applyDarkCurrent()

        if self.interpixCrosstalk:
            #if self.config['verbose'] == 'True': print ("Apply Interpixel Crosstalk")
            print ("\tApply Interpixel Crosstalk")
            self.applyInterpixCrosstalk()

        if self.cosmetics:
            #if self.config['verbose'] == 'True': print ("Add cosmetics")
            print ("\tAdd cosmetics")
            self.addCosmetics()

        if self.nonLinearity:
            #if self.config['verbose'] == 'True': print ("Add non linearity")
            print ("\tAdd non linearity")
            self.applyNonLinearity()

        if self.shutterOpen:
            #if self.config['verbose'] == 'True': print ("Shutter")
            print ("\tShutter")
            self.addReadoutTrails()

        if self.bleeding:
            #if self.config['verbose'] == 'True': print("Apply Saturation")
            print("\tApply Saturation")
            self.applyBleeding()
            
        if self.readoutNoise:
            #if self.config['verbose'] == 'True': print ("Add Readout Noise")
            print ("\tAdd Readout Noise")
            self.applyReadoutNoise()
        """ n'apparait pas, voir tuto, et placer normalement entre ICT et cosmetique
        if self.digiNoise:
            #if self.config['verbose'] == 'True': print ("Add Digitisation Noise")
            print ("\tAdd Digitisation Noise")
            self.applyDigiNoise()
        """
        if self.ADU:
            #if self.config['verbose'] == 'True': print ("electrons2adu")
            print ("\telectrons2adu")
            self.electrons2ADU()

        if self.intscale: 
            #if self.config['verbose'] == 'True': print ("Discretise")
            print ("\tDiscretise")
            self.discretise()

        #if self.config['verbose'] == 'True': print ("Write outputs")
        print ("\tWrite outputs")
        self.writeOutputs()
        
       
        
    def Rampsimulation(self,config_type='file' ):
        #nImage, output_dir,band
        """
        Create a ramp simulated images defined by the configuration file.
        Will do all steps defined in the config file sequentially.

        :return: None
        """
        import copy
        
        for i in range(self.nbGroup*self.nbRead):
            
            # Calcul de la valeur du group et du read pour l'image i
            Group = int(i/self.nbRead)+ 1
            Read = i+1 - self.nbRead*(Group-1)
            
            time= (self.nbRead+self.nbDrop)*(Group-1) + Read
            self.config['exptime'] = time * 1.47528
      
            # Création de l'image vide 
            print ("Read config file and execute ETC")       
            self.configure(config_type)  

            # Set name of output fits file
            if self.acquisition=='ramp':
                self.information['output']='%s/image_%s_M%s_N%s.fits' % (self.output_dir,self.config['filter_band'],Group,Read)
            elif self.acquisition=='CDS' and self.number==1:
                self.information['output']='%s/image_%s_M1_N%s.fits' % (self.output_dir,self.config['filter_band'],Read)
            elif self.acquisition=='CDS' and self.number==2:
                self.information['output']='%s/image_%s_M2_N%s.fits' % (self.output_dir,self.config['filter_band'],Read)
            print ("Building image nb %s: %s:" % (i+1,self.information['output']))                     
                        
            # Set the time between 2 images            
            if Read == 1 and Group != 1 :
                # In this case, if nbDrop!= 0, the time between 2 images is not the exposure time
                self.config['exptime']=1.47528*(self.nbDrop+1)
            else :
                self.config['exptime']=1.47528
                
            # Put the offset or the signal of the n-1 image
            if Read == 1 and Group == 1 : 
                if self.Offset:
                    print ("\tAdd Offset")
                    self.addOffset() 
            else :
                self.image_total+=tabImages[i-1]
 
             
            # Signal ext : ciel + télescope
            if self.addsources:
                #if self.config['verbose'] == 'True': print ("Read objecct list")
                print ("\tGENERATE OBJECTS CATALOG")
                self.generateObjectList()
                self.readObjectlist()
                
                #if self.config['verbose'] == 'True': print ("Generate PSF")
                print ("\tGENERATE PSF")
                self.generatePSF()
                self.readPSFs()
                self.generateFinemaps()
       
                #if self.config['verbose'] == 'True': print ("Add objects")
                print ("\tADD OBJECTS")
                self.addObjects()

            if self.background:
                #if self.config['verbose'] == 'True': print ("Add Sky background")
                print ("\tAdd Sky background")
                self.applySkyBackground()
                
            if self.Vignetting:
                #if self.config['verbose'] == 'True': print ("Add Vignetting")
                print ("\tAdd Vignetting")
                self.applyVignetting()
                    
            if self.cosmicRays:
                #if self.config['verbose'] == 'True': print ("Add cosmic Rays")
                print ("\tAdd cosmic Rays")
                self.addCosmicRays()
              
            # Signal interne : instrument
            if self.persistance:
                #if self.config['verbose'] == 'True': print ("Add dark current")
                print ("\tApply Persistance")
                self.applyPersistance()
                
            if self.instrumIntrinsicNoise:
                #if self.config['verbose'] == 'True': print ("Add dark current")
                print ("\tApply instrument intrinsic noise")
                self.applyInstrumIntrinsicNoise()
            
            if self.darkCurrent:
                #if self.config['verbose'] == 'True': print ("Add dark current")
                print ("\tAdd dark current")
                self.applyDarkCurrent()
            
             
            copie = copy.deepcopy(self.image_total)
            if Read==1 and Group ==1 :
                tabImages=[copie]
            else :
                tabImages.append(copie)
            
            
            # Statistique
            if self.interpixCrosstalk:
                #if self.config['verbose'] == 'True': print ("Apply Interpixel Crosstalk")
                print ("\tApply Interpixel Crosstalk")
                self.applyInterpixCrosstalk()
               
            """ N'apparait pas, voir tuto
            if self.shotNoise:
                #if self.config['verbose'] == 'True': print ("Apply Shot noise")
                print ("\tApply Shot noise")
                self.applyShotNoise()
            """
            
            # Détecteur
            if self.cosmetics:
                #if self.config['verbose'] == 'True': print ("Add cosmetics")
                print ("\tAdd cosmetics")
                self.addCosmetics()
           
            if self.nonLinearity:
                #if self.config['verbose'] == 'True': print ("Add non linearity")
                print ("\tAdd non linearity")
                self.applyNonLinearity()
            
            if self.shutterOpen:
                #if self.config['verbose'] == 'True': print ("Shutter")
                print ("\tShutter")
                self.addReadoutTrails()
            
            if self.bleeding:
                #if self.config['verbose'] == 'True': print("Apply Saturation")
                print("\tApply Saturation")
                self.applyBleeding()
              
            if self.readoutNoise:
                #if self.config['verbose'] == 'True': print ("Add Readout Noise")
                print ("\tAdd Readout Noise")
                self.applyReadoutNoise()
            
            """ n'apparait pas, voir tuto, et placer normalement entre ICT et cosmetique
            if self.digiNoise:
                #if self.config['verbose'] == 'True': print ("Add Digitisation Noise")
                print ("\tAdd Digitisation Noise")
                self.applyDigiNoise()
            """
            if self.ADU:
                #if self.config['verbose'] == 'True': print ("electrons2adu")
                print ("\telectrons2adu")
                self.electrons2ADU()
            
            if self.intscale: 
                #if self.config['verbose'] == 'True': print ("Discretise")
                print ("\tDiscretise")
                self.discretise()
            
            #if self.config['verbose'] == 'True': print ("Write outputs")
            print ("\tWrite outputs")
            self.writeOutputs() 
            

        
    def AcquisitionBruitCDS(self,config_type='file' ):
        """ 
        Créé l'image de bruit CDS
        """
        from math import sqrt
        print("Building CDS image")
        M1_N1=fits.getdata(self.path+'\\images\\AcquisitionBruitCDS\\image_%s_M1_N1.fits' % (self.config['filter_band']))
        M1_N2=fits.getdata(self.path+'\\images\\AcquisitionBruitCDS\\image_%s_M1_N2.fits' % (self.config['filter_band']))
        image1 = np.float32(M1_N2) - np.float32(M1_N1)
        M2_N1=fits.getdata(self.path+'\\images\\AcquisitionBruitCDS\\image_%s_M2_N1.fits' % (self.config['filter_band']))
        M2_N2=fits.getdata(self.path+'\\images\\AcquisitionBruitCDS\\image_%s_M2_N2.fits' % (self.config['filter_band']))
        image2 = np.float32(M2_N2) - np.float32(M2_N1)
        
        gain_conv=fits.getdata(self.path+'/data/GainMap/'+self.information['GainMapFile'])
    
        CDSNoise = (image2 - image1)/sqrt(2)*gain_conv
        hdu=fits.PrimaryHDU()
        hdu.data=CDSNoise
        hdu.writeto(self.path+'\\images\\AcquisitionBruitCDS\\CDSNoise.fits',overwrite=True)
        
        print("end")


if __name__ == '__main__':

    #run the simulator
    IS = ImageSimulator()
    IS.simulate()



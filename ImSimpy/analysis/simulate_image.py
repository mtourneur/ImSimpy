import os
import sys
import subprocess as sp
import numpy as np
from astropy.io import ascii
from astropy.table import Table

""" This program will generate astronomical image for the GFT with Skymaker, extract the source
    with sextractor and study the completness.
    It assumes that the PSF is already done, with the format 512x512 and a pixel size of 3 microns
    Requires Skymaker and SExtractor to be installed
    Python3.5
"""

def delete_nearby_stars(cat,limit):
   """ delete sources within the limit distance """
   # copy the input catalogue
   new_cat=cat.copy()
   rows_to_delete=[]
   i1=0
   for star in cat:
       i2=0
       for star2 in cat:
           if i1 not in rows_to_delete and i2!=i1:

               dist=np.sqrt((star2['COORD_XPIXEL']-star['COORD_XPIXEL'])**2+(star2['COORD_YPIXEL']-star['COORD_YPIXEL'])**2)

               if dist < limit:
                   rows_to_delete.append(i2)
           i2+=1
       i1+=1
   # Delete duplicates in the list, and sort it
   rows_to_delete=sorted(set(rows_to_delete))
   new_cat.remove_rows(rows_to_delete)
   print ("Number of closed stars deleted: {0:d}\n".format(len(rows_to_delete)))
   return new_cat



def stars_generator(N_stars,xlims,ylims,maglims,output_fname,close_stars,min_dist):
   """ Generate randomly the catalogue of stars  """
   import numpy.random as rdm
 
   # Generate random set of parameters from uniform distribution within intervals"
   x1,x2=xlims
   y1,y2=ylims
   m1,m2=maglims

   x = rdm.uniform(x1,x2,int(N_stars))
   y = rdm.uniform(y1,y2,int(N_stars))
   mag = rdm.uniform(m1,m2,int(N_stars))

   # Add the object code for having the same format as a stuff catalogue
   obj_code=np.ones(int(N_stars))*100

   #Write position in text file
   np.savetxt('stars_cat.list', list(zip(obj_code,x,y,mag)),fmt='%d %.2f %.2f %.2f')

   # If file already exists delete the catalogue of sources
   try:
       os.remove('%s.list' % output_fname)
   except OSError:
       pass

   #Add the Headers for the stuff catalogue in order to read it as a sextractor file by astropy.ascii
   sp.run('cat stuff_headers_stars.txt stars_cat.list >> %s.list' % output_fname,shell=True)

   # Delete the 'stars_cat.list' which is now useless
   #os.remove('stars_cat.list')

   if close_stars == False:
       print ("\nDeleting stars closer than {0:d} pixels...\n".format(min_dist))
       stars_cat=ascii.read('%s.list' % output_fname)
       new_stars_cat=delete_nearby_stars(stars_cat,min_dist)
       ascii.write(new_stars_cat,'stars_cat_notclosed.list')

       # Delete First line 
       sp.run("sed -i '1d' stars_cat_notclosed.list",shell=True)

       try:
           os.remove('%s.list' % output_fname)
       except OSError:
           pass

       #Add the Headers for the stuff catalogue in order to read it as a sextractor file by astropy.ascii
       sp.run('cat stuff_headers_stars.txt stars_cat_notclosed.list >> %s.list' % output_fname,shell=True)

       #os.remove('stars_cat_notclosed.list')

def mag_distrib(sources_cat,mag_lims):
   """ Substitute the magnitudes from the input sources catalog """
   import numpy.random as rdm

   stars_cat=ascii.read('sources_catalogs/stars_cat_notclosed.list')
   m1,m2=mag_lims
  
   mag = rdm.uniform(m1,m2,len(stars_cat))
   stars_cat['col4']=[round(x,2) for x in mag] 
   ascii.write(stars_cat,'sources_catalogs/temp.list',comment=False)

   # Delete first line of the file. It corresponds to the headers of the astropy Table
   sp.run("sed -i '1d' sources_catalogs/temp.list",shell=True) 

   #Add the Headers for the stuff catalogue in order to read it as a sextractor file by astropy.ascii
   sp.run('cat sources_catalogs/stuff_headers_stars.txt sources_catalogs/temp.list > %s.list' % sources_cat,shell=True)
   
   try:
       os.remove('sources_catalogs/temp.list')
   except OSError:
       pass
   

def create_image(sources_cat,skymaker_pars):
   """ Use skymaker to simulate images """

   
   ImageSize,image_fname,ImageType,gain,well_cap,satur_level,RN,EXPTIME,ZP,pixelSize,psf_fname,psf_oversamp,psf_map,DM1,DM2,wvl,SB,nthreads  = skymaker_pars

   sp.run ('sky %s.list -c gft.sky -IMAGE_NAME %s.fits -IMAGE_SIZE %d -IMAGE_TYPE %s -GAIN %f -WELL_CAPACITY %f -SATUR_LEVEL %d -READOUT_NOISE %f -EXPOSURE_TIME %f -MAG_ZEROPOINT %f -PIXEL_SIZE %f -PSF_TYPE FILE -PSF_NAME %s.fits -PSF_OVERSAMP %f -PSF_MAPSIZE %f -M1_DIAMETER %f -M2_DIAMETER %f -WAVELENGTH %f -BACK_MAG %f -STARCOUNT_ZP 0 -VERBOSE_TYPE FULL -NTHREADS %d' % (sources_cat,image_fname,ImageSize[0],ImageType,gain,well_cap,satur_level,RN,EXPTIME,ZP,pixelSize,psf_fname,psf_oversamp,psf_map,DM1,DM2,wvl,SB,nthreads) ,shell=True)
  

def sources_extraction(image,sextractor_pars):
   """ Extract sources from the generated image using sextractor """

   cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP, gain, pixelScale,seeing,back_type,back_value,back_size,backphoto_type,backphoto_thick,back_filterthresh,checkimage_type,checkimage_name= sextractor_pars
   sp.run('sex %s.fits -c gft.sex -CATALOG_NAME %s.cat -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME gft.param -DETECT_TYPE CCD -DETECT_MINAREA %d -DETECT_THRESH %d -ANALYSIS_THRESH %d -PHOT_APERTURES %d -SATUR_LEVEL %d -MAG_ZEROPOINT %f -GAIN %f -PIXEL_SCALE %f -SEEING_FWHM %f -BACK_TYPE %s -BACK_VALUE %f -BACK_SIZE %d -BACKPHOTO_TYPE %s -BACKPHOTO_THICK %d -BACK_FILTTHRESH %f -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s.fits ' % (image,cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP, gain, pixelScale,seeing,back_type,back_value,back_size,backphoto_type,backphoto_thick,back_filterthresh,checkimage_type,checkimage_name),shell=True)

   

def completeness(input_sources_cat,detected_sources_cat,output_fname,Mag_lim,pix_radius):
   """ Finds out how many stars were detected """
   from astropy.table import Column

   #Load catalogues in table
   input_cat=ascii.read('%s.list' % input_sources_cat)
   detected_cat=ascii.read('%s.cat' % detected_sources_cat) 

   print ('Number of sources in stuff catalog above the mag lim of %.2f: %d' % (Mag_lim,len(input_cat[input_cat['MAG_APP']<Mag_lim])))
   print ('Number of sources detected: %d \n' % len(detected_cat))

   #Pixel radius
   pixradius=pix_radius

   nb=0
   i=0
   det=np.zeros(len(input_cat))
   x_det_list=np.zeros(len(input_cat))
   y_det_list=np.zeros(len(input_cat))
   mag_sex=np.zeros(len(input_cat))

   col_det=Column(name='detected',data=det)
   x_det_coord=Column(name='x_coord_det',data=x_det_list)
   y_det_coord=Column(name='y_coord_det',data=y_det_list)
   mag_det=Column(name='mag_det',data=mag_sex)
   input_cat.add_columns([col_det,x_det_coord,y_det_coord,mag_det])

   col_det_sex=Column(name='detected',data=np.zeros(len(detected_cat)))
   detected_cat.add_columns([col_det_sex])


   for x1, y1 in zip (detected_cat['XPEAK_IMAGE'], detected_cat['YPEAK_IMAGE']):
       #print ('object n. {0:d} at position: {1:.2f}-{2:.2f} \n'.format(nb,x1,y1))
       min_dist=1e40
       j=0
       x_det=-1;y_det=-1;
       for x2,y2,mag in zip(input_cat['COORD_XPIXEL'],input_cat['COORD_YPIXEL'],input_cat['MAG_APP']):
           if detected_cat['detected'][i]==0 and x1 >= int(x2)-pixradius and x1 <= int(x2)+pixradius and y1 >= int(y2)-pixradius and y1 <= int(y2)+pixradius:
               #Test the minimum distance
               dist=(x2-x1)**2+(y2-y1)**2
               if dist < min_dist and detected_cat['MAG_AUTO'][i] > 0.9*mag and detected_cat['MAG_AUTO'][i] < 1.1*mag:
                  min_dist=dist
                  x_det=x1
                  y_det=y1
                  mag_det=detected_cat['MAG_AUTO'][i]
                  index=j
           j+=1
       if min_dist<1e40:
           nb+=1
           detected_cat['detected'][i]=1
           #print ('Matched sources n. {0:d} at position: {1:.2f}-{2:.2f} \n'.format(i,x_det,y_det))
           input_cat['detected'][index]=1
           input_cat['x_coord_det'][index]=x_det
           input_cat['y_coord_det'][index]=y_det
           input_cat['mag_det'][index]=mag_det
       else:
           detected_cat['detected'][i]=-1
           #print ('Matched sources n. {0:d} at position: {1:.2f}-{2:.2f} \n'.format(i,x_det,y_det))

       i+=1


   """
   for x1,y1 in zip(input_cat['COORD_XPIXEL'],input_cat['COORD_YPIXEL']):
       nb+=1
       #print ('object n. {0:d} at position: {1:.2f}-{2:.2f} \n'.format(nb,x1,y1))
       min_dist=1e40
       x_det=-1;y_det=-1;
       j=0
       for x2, y2 in zip (detected_cat['XPEAK_IMAGE'], detected_cat['YPEAK_IMAGE']):
           if detected_cat['detected'][j]==0 and x2 >= int(x1)-pixradius and x2 <= int(x1)+pixradius and y2 >= int(y1)-pixradius and y2 <= int(y1)+pixradius:
               #Test the minimum distance
               dist=(x2-x1)**2+(y2-y1)**2
               if dist < min_dist:
                  min_dist=dist
                  x_det=x2
                  y_det=y2
                  mag_det=detected_cat['MAG_AUTO'][j]
                  index=j
           j+=1
   
       if min_dist<1e40:
            i+=1
            detected_cat['detected'][index]=1
            #print ('Matched sources n. {0:d} at position: {1:.2f}-{2:.2f} \n'.format(i,x_det,y_det))
            input_cat['detected'][nb-1]=1
            input_cat['x_coord_det'][nb-1]=x_det
            input_cat['y_coord_det'][nb-1]=y_det
            input_cat['mag_det'][nb-1]=mag_det
   """
   #Cross match catalog
   print ('Number of sources matched in both catalogs: %d' % nb)

   #Write output file
   ascii.write(input_cat,'%s.txt' % output_fname)



   x_false_list=detected_cat['XPEAK_IMAGE'][detected_cat['detected']==-1]
   y_false_list=detected_cat['YPEAK_IMAGE'][detected_cat['detected']==-1]
   mag_sex=detected_cat['MAG_AUTO'][detected_cat['detected']==-1]

   #x_det_coord=Column(name='x_coord',data=x_det_list)
   #y_det_coord=Column(name='y_coord',data=y_det_list)
   #mag_det=Column(name='mag_det',data=mag_sex)
   false_det_cat=Table([x_false_list,y_false_list,mag_sex],names=('x_coord','y_coord','mag_det'))


   #Write false detections in a separated file
   ascii.write(false_det_cat,'%s_false_detections.txt' % output_fname)


def plot_completeness(cat_name,name_plot,mag_lims,binning_mag,second_cat='no'):
   """ Plot the number of detected sources vs MAgnitude """

   cat=ascii.read('%s.txt' % cat_name)
   mag_bins=np.arange(mag_lims[0],mag_lims[1],binning_mag)

   mask=cat['detected']==1
   mag_binned_tot=np.digitize(cat['MAG_APP'],mag_bins,right=True)
   mag_binned_det=np.digitize(cat[mask]['MAG_APP'],mag_bins,right=True)

   nb_mag=np.array([  len(np.where(mag_binned_tot==i)[0]) for i in range(1,len(mag_bins)) ])
   nb_mag_det = np.array([  len(np.where(mag_binned_det==i)[0]) for i in range(1,len(mag_bins)) ])
   #mag_tot= np.array([stuff_cat['MAG_APP'][mag_binned_tot == i].mean() for i in range(1, len(mag_bins))])
    #mag_det= np.array([stuff_cat[mask]['MAG_APP'][mag_binned_det == i].mean() for i in range(1, len(mag_bins))])
   print (nb_mag)
   print (nb_mag_det)

   #Write completeness result in text file
   np.savetxt('completeness_table.csv', list(zip(mag_bins,nb_mag,nb_mag_det)),fmt='%.2f %d %d')


   mag_bin_plot=(mag_bins[:-1]+mag_bins[1:])/2

   import matplotlib.pyplot as plt

   # the histogram of the input sources
   n, bins, patches = plt.hist(cat['MAG_APP'], mag_bins, normed=0, facecolor='green', alpha=0.75)
   plt.xlabel('Magnitude')
   plt.ylabel('Nb of sources')
   plt.xlim([mag_bins[0],mag_bins[-1]])
   plt.savefig('plots/hist_sources.png')
   #plt.show()

   plt.clf()
   plt.plot(mag_bin_plot,nb_mag_det/nb_mag)
   plt.xlabel('Magnitude AB')
   plt.ylabel('Efficiency')
   plt.grid(True)
   plt.savefig('plots/completeness_%s.png' % name_plot)
   plt.show()


   if second_cat != 'no':
       cat2=ascii.read('%s.txt' % second_cat)
       mag_bins2=np.arange(mag_lims[0],mag_lims[1],binning_mag)

       mask2=cat2['detected']==1
       mag_binned_tot2=np.digitize(cat2['MAG_APP'],mag_bins2,right=True)
       mag_binned_det2=np.digitize(cat2[mask2]['MAG_APP'],mag_bins2,right=True)

       nb_mag2=np.array([  len(np.where(mag_binned_tot2==i)[0]) for i in range(1,len(mag_bins2)) ])
       nb_mag_det2 = np.array([  len(np.where(mag_binned_det2==i)[0]) for i in range(1,len(mag_bins2)) ])

       mag_bin_plot2=(mag_bins2[:-1]+mag_bins2[1:])/2
       #print (mag_bin_plot)
       #plt.plot(mag_bin_plot,nb_mag_det/nb_mag,label='seeing=0.7"',color='red')
       #plt.plot(mag_bin_plot2,nb_mag_det2/nb_mag2,label='seeing=0.1"',color='green')
       plt.plot(mag_bin_plot,nb_mag_det/nb_mag,label='5.9',color='red')
       plt.plot(mag_bin_plot2,nb_mag_det2/nb_mag2,label='5',color='green')
       plt.xlabel('Magnitude AB')
       plt.ylabel('Efficiency')
       #plt.yscale('log')
       #plt.xscale('log')
       plt.grid(True)
       plt.legend()
       plt.savefig('plots/completeness_comp.png')
       plt.show()


def plot_photometric_accuracy(cat_name,name_plot):
   """ Plot the magnitude returned by SExtractor vs the input magnitude """

   cat=ascii.read('%s.txt' % cat_name)
   mask = cat['detected']==1
   
   import matplotlib.pyplot as plt
   x=np.linspace(min(cat[mask]['MAG_APP']),max(cat[mask]['MAG_APP']),1000)
   plt.clf()
   plt.xlim(int(min(cat[mask]['MAG_APP']))-1,np.ceil(max(cat[mask]['MAG_APP']))+1)
   plt.ylim(int(min(cat[mask]['MAG_APP']))-1,np.ceil(max(cat[mask]['MAG_APP']))+1)
   plt.xticks(np.arange(int(min(cat[mask]['MAG_APP']))-1,int(np.ceil(max(cat[mask]['MAG_APP'])))+1))
   
   plt.yticks(np.arange(int(min(cat[mask]['MAG_APP']))-1,int(np.ceil(max(cat[mask]['MAG_APP'])))+1))
   plt.scatter(cat[mask]['MAG_APP'],cat[mask]['mag_det'])
   plt.plot(x,x,color='red',ls='--')
   plt.xlabel('Input Magnitude (AB)')
   plt.ylabel('Extracted magnitude (AB)')
   plt.grid(True)
   plt.savefig('plots/phot_accuracy_%s.png' % name_plot)
   plt.show()


def create_ds9_regions(cat_name,output_fname,color,radius,Type):
   """ Create a ds9 region, in our case it just consists in circling the sources """
   cat=ascii.read('%s.txt' % cat_name)
   with open('%s.reg' % output_fname, 'w') as f:
       f.write('global color=%s dashlist=8 3 width=1\n' % color)
       f.write('image\n')
       if Type == 'input':
           for x,y,mag in zip(cat['COORD_XPIXEL'],cat['COORD_YPIXEL'],cat['MAG_APP']):
               f.write('circle(%d,%d,%d)  # text = {%.2f}\n' %(x,y,radius,mag))
       elif Type=='detected':
           for x,y,mag in zip(cat['x_coord_det'],cat['y_coord_det'],cat['mag_det']):
               f.write('circle(%d,%d,%d)  # text = {                   %.2f}\n' %(x,y,radius,mag))
       elif Type=='false_detection':
           for x,y,mag in zip(cat['x_coord'],cat['y_coord'],cat['mag_det']):
               f.write('circle(%d,%d,%d)  # text = {                   %.2f}\n' %(x,y,radius,mag))


if __name__ == '__main__':
 

   generate_sources=0
   change_mag=1
   generate_image=1
   do_source_extraction=1
   do_completeness=1
   display=0
  

   results_dir='sources_catalogs/'


   # Define parameters for sources generation
   zN_stars=5000
   Xlims = [100,3900]
   Ylims = [100,3900]
   input_cat_name=results_dir+'input_sources'
   allow_close_stars=False
   min_distance_between_2stars=100

   # Random generator seed
   seed = 2021
   np.random.seed(seed)

   # Instrument parameters
   #gain = 2.0#5.34
   #well_cap = 350000.0
   satur_level = 65535
   #RN=8     # 8    12
   #ZP=23.79  #21.54   # 23.29
   #pixelScale=0.381      # in microns
   DM1=1.3
   DM2=0.47
   #SB= 18.7             # in mag/arcsec2
   EXPTIME=600         # in seconds
   seeing=0.7            # FWHM in arcsecond  between 0.025 and 5.5"

   band='H'
   if band == 'g':
       gain=5.34
       ZP=24.60 -2.5*np.log10(gain)
       SB=20.55
       pixelScale=0.381
       psf_oversamp=10.0267
       wvl=0.480
       well_cap=350000
       RN=8
       ImageSize=[4096,4096]
       mag_lims=[20,26]
   elif band == 'z':
       gain=5.34
       ZP=23.79-2.5*np.log10(gain)
       SB=18.7
       pixelScale=0.381
       psf_oversamp=5.5473
       wvl=0.867
       well_cap=350000
       RN=8
       ImageSize=[4096,4096]
       mag_lims=[18,24]
    
   elif band == 'H':
       gain=2.0
       ZP=23.94 -2.5*np.log10(gain)
       SB=15.4
       pixelScale=0.762
       psf_oversamp=5.8747
       wvl=1.6368
       well_cap=80000
       RN=24
       EXPTIME=30
       ImageSize=[4096,4096]
       mag_lims=[16,22]
       gain*=20
       well_cap*=20


   # Skymaker parameters
   ImageSize=[4096,4096]
   ImageType='SKY'      # 'SKY' or 'SKY_NONOISE'
   image_dir='images/'
   image_fname='zemax_zband1'
   psf_name='psf/skymaker/psf_zemax_H_edge_moffat_s10_128_3064mic_50'
   #psf_oversamp=5.5473
   psf_map=256
   #wvl=0.867           # in microns
   nthreads=1
   
   image_fname=image_fname+'_'+ImageType.lower()
 
   skymaker_pars=ImageSize,image_dir+image_fname,ImageType,gain,well_cap,satur_level,RN,EXPTIME,ZP,pixelScale,psf_name,psf_oversamp,psf_map,DM1,DM2,wvl,SB,nthreads
   #skymaker_pars=ImageSize,image_fname,ImageType,gain,well_cap,RN,EXPTIME,ZP,pixelScale,psf_name,psf_oversamp,psf_map,DM1,DM2,wvl,SB,nthreads


   #Sextractor parameters
   detected_cat_name=results_dir+'zemax_zband_sex'
   detect_minarea = 1
   detect_thresh = 5
   analysis_thresh = 1 
   phot_aperture = 10
   back_type = 'AUTO'
   back_value = 0.0
   back_size = 256
   ZP_sextractor = ZP + 2.5*np.log10(EXPTIME)
   backphoto_type = 'GLOBAL'
   backphoto_thick = 24
   backphoto_filterthresh = 0.0
   checkimage_type = 'NONE'
   checkimage_name = 'all_sources'
    
   checkimage_name = image_fname +'_'+ checkimage_name

   sextractor_pars = detected_cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP_sextractor, gain, pixelScale,seeing,back_type,back_value,back_size,backphoto_type,backphoto_thick,backphoto_filterthresh,checkimage_type,checkimage_name

   # Parameters for completeness comparison
   pix_radius = 5      # search for sources in: pix_input_source +/- pix_radius 
   Mag_lim = 22.64     # From etc
   binning_mag=0.1     # bin of 0.5 between the mag_lims


   # Parameters for ds9 regions
   plots_dir='plots/'
   cat_crossmatch=results_dir+'crossmatch'
   radius=10      # diameter of the circle in pixels
   factor=1.3     # radius * factor: radius of the dectected source in pixels
   color_input='red'
   color_detected='green'
   color_false='blue'
   output_dir='images/region_ds9/'
   output_name_in='input'
   output_name_det='detected'
   output_name_false='false_detection'
   false_detection_name='false_detections'

   output_name_in=output_dir+output_name_in+'_'+image_fname
   output_name_det=output_dir+output_name_det+'_'+image_fname
   output_name_false=output_dir+output_name_false+'_'+image_fname
   cat_crossmatch=cat_crossmatch+'_'+image_fname
   cat_false_detections=cat_crossmatch+'_'+false_detection_name

   if generate_sources:
       # Generate the stars catalogue
       stars_generator(N_stars,Xlims,Ylims,mag_lims,input_cat_name,allow_close_stars,min_distance_between_2stars)

   if change_mag:
       mag_distrib(input_cat_name,mag_lims)

   if generate_image:
       # Create image
       create_image(input_cat_name,skymaker_pars)

   if do_source_extraction:
       # extract sources of the simulated image with sextractor
       sources_extraction(image_dir+image_fname,sextractor_pars)

   if do_completeness:
       # Compares how many stars were detected with the input catalogue
       completeness(input_cat_name,detected_cat_name,cat_crossmatch,Mag_lim,pix_radius) 
       # Plot the nb of detected sources vs magnitude
       plot_completeness(cat_crossmatch,image_fname,mag_lims,binning_mag,second_cat=results_dir+'crossmatch_zemax_zband5_sky')
       # Plot the photometric accuracy
       plot_photometric_accuracy(cat_crossmatch,image_fname)

       # create ds9 regions to circle input sources
       create_ds9_regions(cat_crossmatch,output_name_in,color_input,radius,'input')
       # create ds9 regions to circle input sources
       create_ds9_regions(cat_crossmatch,output_name_det,color_detected,radius*1.3,'detected') 
       # create ds9 regions to circle false detections
       create_ds9_regions(cat_false_detections,output_name_false,color_false,radius*factor*1.5,'false_detection') 

 
   if display:
       # Show the simulated image with ds9 regions indicating input and detected sources
       # Uses ds9
       sp.run('ds9 %s.fits -region %s.reg -region %s.reg -region %s.reg' % (image_dir+image_fname,output_name_in,output_name_det,output_name_false),shell=True) 

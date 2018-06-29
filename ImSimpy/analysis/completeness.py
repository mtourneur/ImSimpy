import os
import sys
import subprocess as sp
import numpy as np
from astropy.io import ascii
from astropy.table import Table,Column

def sources_extraction(image,sextractor_pars):
   """ Extract sources from the generated image using sextractor """

   cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP, gain, pixelScale,seeing,back_type,back_value,back_size,backphoto_type,backphoto_thick,back_filterthresh,checkimage_type,checkimage_name= sextractor_pars
   sp.run('sex %s.fits -c gft.sex -CATALOG_NAME %s.cat -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME gft.param -DETECT_TYPE CCD -DETECT_MINAREA %d -DETECT_THRESH %d -ANALYSIS_THRESH %d -PHOT_APERTURES %d -SATUR_LEVEL %d -MAG_ZEROPOINT %f -GAIN %f -PIXEL_SCALE %f -SEEING_FWHM %f -BACK_TYPE %s -BACK_VALUE %f -BACK_SIZE %d -BACKPHOTO_TYPE %s -BACKPHOTO_THICK %d -BACK_FILTTHRESH %f -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s.fits ' % (image,cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP, gain, pixelScale,seeing,back_type,back_value,back_size,backphoto_type,backphoto_thick,back_filterthresh,checkimage_type,checkimage_name),shell=True)

   

def completeness(input_sources_cat,detected_sources_cat,output_fname,cat_falsedet,Mag_lim,pix_radius):
   """ Finds out how many stars were detected """

   #Load catalogues in table
   input_cat=ascii.read('%s.txt' % input_sources_cat)
   detected_cat=ascii.read('%s.cat' % detected_sources_cat) 
   #print (input_cat)
   #print (detected_cat)
   print ('Number of sources in stuff catalog below the mag lim of %.2f: %d' % (Mag_lim,len(input_cat[input_cat['MAG']<Mag_lim])))
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
       for x2,y2,mag in zip(input_cat['COORD_XPIXEL'],input_cat['COORD_YPIXEL'],input_cat['MAG']):
           if detected_cat['detected'][i]==0 and x1 >= int(x2)-pixradius and x1 <= int(x2)+pixradius and y1 >= int(y2)-pixradius and y1 <= int(y2)+pixradius:
               #Test the minimum distance
               dist=(x2-x1)**2+(y2-y1)**2
               if dist < min_dist:# and detected_cat['MAG_AUTO'][i] > 0.9*mag and detected_cat['MAG_AUTO'][i] < 1.1*mag:
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
   for x1,y1 in zip(input_cat['COORD_YPIXEL'],input_cat['COORD_XPIXEL']):
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
   ascii.write(false_det_cat,'%s.txt' % cat_falsedet)


def plot_completeness(cat_name,output_name,name_plot,mag_lims,binning_mag,plot,second_cat='no'):
   """ Plot the number of detected sources vs MAgnitude """

   cat=ascii.read('%s.txt' % cat_name)
   mag_bins=np.arange(mag_lims[0],mag_lims[1],binning_mag)

   mask=cat['detected']==1
   mag_binned_tot=np.digitize(cat['MAG'],mag_bins,right=True)
   mag_binned_det=np.digitize(cat[mask]['MAG'],mag_bins,right=True)

   nb_mag=np.array([  len(np.where(mag_binned_tot==i)[0]) for i in range(1,len(mag_bins)) ])
   nb_mag_det = np.array([  len(np.where(mag_binned_det==i)[0]) for i in range(1,len(mag_bins)) ])
   #mag_tot= np.array([stuff_cat['MAG'][mag_binned_tot == i].mean() for i in range(1, len(mag_bins))])
    #mag_det= np.array([stuff_cat[mask]['MAG'][mag_binned_det == i].mean() for i in range(1, len(mag_bins))])
   print (nb_mag)
   print (nb_mag_det)

   #Write completeness result in text file
   np.savetxt('%s.txt' % output_name, list(zip(mag_bins,nb_mag,nb_mag_det)),fmt='%.2f %d %d')


   mag_bin_plot=(mag_bins[:-1]+mag_bins[1:])/2

   import matplotlib.pyplot as plt

   # the histogram of the input sources
   n, bins, patches = plt.hist(cat['MAG'], mag_bins, normed=0, facecolor='green', alpha=0.75)
   plt.xlabel('Magnitude')
   plt.ylabel('Nb of sources')
   plt.xlim([mag_bins[0],mag_bins[-1]])
   plt.savefig('results/plots/hist_sources.png')
   #plt.show()

   plt.clf()
   plt.plot(mag_bin_plot,nb_mag_det/nb_mag)
   plt.xlabel('Magnitude AB')
   plt.ylabel('Efficiency')
   plt.grid(True)
   plt.savefig('%s.png' % output_name)
   if plot: plt.show()


   if second_cat != 'no':
       cat2=ascii.read('%s.txt' % second_cat)
       mag_bins2=np.arange(mag_lims[0],mag_lims[1],binning_mag)

       mask2=cat2['detected']==1
       mag_binned_tot2=np.digitize(cat2['MAG'],mag_bins2,right=True)
       mag_binned_det2=np.digitize(cat2[mask2]['MAG'],mag_bins2,right=True)

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
       plt.savefig('results/plots/completeness_comp.png')
       if plot: plt.show()


def plot_photometric_accuracy(cat_name,name_plot,plot):
   """ Plot the magnitude returned by SExtractor vs the input magnitude """

   cat=ascii.read('%s.txt' % cat_name)
   mask = cat['detected']==1
   
   import matplotlib.pyplot as plt
   x=np.linspace(min(cat[mask]['MAG']),max(cat[mask]['MAG']),1000)
   plt.clf()
   plt.xlim(int(min(cat[mask]['MAG']))-1,np.ceil(max(cat[mask]['MAG']))+1)
   plt.ylim(int(min(cat[mask]['MAG']))-1,np.ceil(max(cat[mask]['MAG']))+1)
   plt.xticks(np.arange(int(min(cat[mask]['MAG']))-1,int(np.ceil(max(cat[mask]['MAG'])))+1))
   
   plt.yticks(np.arange(int(min(cat[mask]['MAG']))-1,int(np.ceil(max(cat[mask]['MAG'])))+1))
   plt.scatter(cat[mask]['MAG'],cat[mask]['mag_det'])
   plt.plot(x,x,color='red',ls='--')
   plt.xlabel('Input Magnitude (AB)')
   plt.ylabel('Extracted magnitude (AB)')
   plt.grid(True)
   plt.savefig('results/plots/%s.png' % name_plot)
   if plot: plt.show()


def create_ds9_regions(cat_name,output_fname,color,radius,Type):
   """ Create a ds9 region, in our case it just consists in circling the sources """
   cat=ascii.read('%s.txt' % cat_name)
   with open('%s.reg' % output_fname, 'w') as f:
       f.write('global color=%s dashlist=8 3 width=1\n' % color)
       f.write('image\n')
       if Type == 'input':
           for x,y,mag in zip(cat['COORD_XPIXEL'],cat['COORD_YPIXEL'],cat['MAG']):
               f.write('circle(%d,%d,%d)  # text = {%.2f}\n' %(x,y,radius,mag))
       elif Type=='detected':
           for x,y,mag in zip(cat['x_coord_det'],cat['y_coord_det'],cat['mag_det']):
               f.write('circle(%d,%d,%d)  # text = {                   %.2f}\n' %(x,y,radius,mag))
       elif Type=='false_detection':
           for x,y,mag in zip(cat['x_coord'],cat['y_coord'],cat['mag_det']):
               f.write('circle(%d,%d,%d)  # text = {                   %.2f}\n' %(x,y,radius,mag))


if __name__ == '__main__':
   import time
 
   do_source_extraction=1
   do_completeness=1
   display=1
 
   #Sextractor parameters
   image_dir='../images/'#'../images/psf_degrad_30s/'
   detect_minarea = 1
   detect_thresh = 4
   analysis_thresh = 1
   phot_aperture = 10
   back_type = 'AUTO'
   back_value = 0.0
   back_size = 256
   backphoto_type = 'GLOBAL'
   backphoto_thick = 24
   backphoto_filterthresh = 0.0
   checkimage_type = 'NONE'
   checkimage_name = 'all_sources'

   # Instrument parameters
   satur_level = 65535
   
   # Parameters for completeness comparison
   pix_radius = 5      # search for sources in: pix_input_source +/- pix_radius 
   binning_mag=0.1     # bin of 0.5 between the mag_lims

   plot=False
   # Parameters for ds9 regions
   output_dir='test/'
   radius=10      # diameter of the circle in pixels
   factor=1.3     # radius * factor: radius of the dectected source in pixels
   color_input='red'
   color_detected='green'
   color_false='blue'
   ds9Reg_dir='region_ds9/'


   bands=['i']
   locs=['ideal']
   seeings=['07']#,'07','10']            
   zoom=[1]#[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,3.0,4.0]  #[None]
    
   for band in bands:
       input_cat_name='../data/catalog/input_sources_%s' % band +'_30s'
       for loc in locs:
           for seeing in seeings:
               if zoom[0] != None:
                   for z in zoom:
                       image_fname='science_%s_%s_s%s_z%s' % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       detected_cat_name='results/sex_detected_cat_%s_%s_s%s_z%s' % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       cat_crossmatch='results/crossmatch_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       cat_false_detections='results/false_detections_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       completeness_file='results/completeness_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       phot_accuracy_plot='phot_accuracy_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       ds9Regions_input=ds9Reg_dir+'%s_%s_s%s_z%s_input'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       ds9Regions_det=ds9Reg_dir+'%s_%s_s%s_z%s_det'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       ds9Regions_falsedet=ds9Reg_dir+'%s_%s_s%s_z%s_falsedet'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                       if band == 'g':
                           gain=1.5
                           ZP=24.3 -2.5*np.log10(gain)
                           pixelScale=0.381
                           well_cap=350000
                           EXPTIME=600
                           mag_lims=[18,24]
                           lim_mag = 22.15     # From etc
                       elif band == 'r':
                           gain=1.5
                           ZP=24.24-2.5*np.log10(gain)
                           pixelScale=0.381
                           well_cap=350000
                           EXPTIME=30
                           mag_lims=[16,22]
                           lim_mag = 22.15     # From etc
                       elif band == 'i':
                           gain=1.5
                           ZP=23.94-2.5*np.log10(gain)
                           pixelScale=0.381
                           well_cap=350000
                           EXPTIME=30
                           mag_lims=[16,22]
                           lim_mag = 22.15     # From etc

                       elif band == 'z':
                           gain=1.5
                           ZP=23.24-2.5*np.log10(gain)
                           pixelScale=0.381
                           well_cap=350000
                           EXPTIME=30
                           mag_lims=[16,22]
                           lim_mag = 22.15     # From etc
                       elif band == 'H':
                           gain=1.5 
                           ZP=23.94 -2.5*np.log10(gain)
                           pixelScale=0.762
                           well_cap=80000
                           EXPTIME=30*20
                           mag_lims=[16,22]
                           lim_mag = 22.15     # From etc

                       ZP_sextractor = ZP + 2.5*np.log10(EXPTIME)
                       checkimage_name = image_fname +'_'+ checkimage_name

                       sextractor_pars = detected_cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP_sextractor, gain, pixelScale,float(seeing)/10,back_type,back_value,back_size,backphoto_type,backphoto_thick,backphoto_filterthresh,checkimage_type,checkimage_name

                       time.sleep(1)
                       if do_source_extraction:
                           # extract sources of the simulated image with sextractor
                           sources_extraction(image_dir+image_fname,sextractor_pars)

                       if do_completeness:
                           # Compares how many stars were detected with the input catalogue
                           completeness(input_cat_name,detected_cat_name,cat_crossmatch,cat_false_detections,lim_mag,pix_radius) 
                           # Plot the nb of detected sources vs magnitude
                           plot_completeness(cat_crossmatch,completeness_file,image_fname,mag_lims,binning_mag,plot)
                           # Plot the photometric accuracy
                           plot_photometric_accuracy(cat_crossmatch,phot_accuracy_plot,plot)

                           # create ds9 regions to circle input sources
                           create_ds9_regions(cat_crossmatch,ds9Regions_input,color_input,radius,'input')
                           # create ds9 regions to circle input sources
                           create_ds9_regions(cat_crossmatch,ds9Regions_det,color_detected,radius*1.3,'detected') 
                           # create ds9 regions to circle false detections
                           create_ds9_regions(cat_false_detections,ds9Regions_falsedet,color_false,radius*factor*1.5,'false_detection') 

 
                       if display:
                            # Show the simulated image with ds9 regions indicating input and detected sources
                            # Uses ds9
                            sp.run('ds9 %s.fits -region %s.reg -region %s.reg -region %s.reg' % (image_dir+image_fname,ds9Regions_input,ds9Regions_det,ds9Regions_falsedet),shell=True)
                   

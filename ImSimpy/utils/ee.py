import numpy as np
from photutils import CircularAperture,aperture_photometry
from astropy.io import fits

import matplotlib.pyplot as plt
import sys

psf_file=[]

for arg in sys.argv:
   psf_file.append(arg)

#print (psf_file[-1].split(), psf_file[-1].split()[0])

band='z'
if band =='g':
   pixel_size_real=15
   pixel_size=0.187   # in microns
   radii = np.arange(1,850)
elif band =='z':
   pixel_size_real=15
   pixel_size=0.338   # in microns
   radii = np.arange(1,800)
elif band =='H':
   pixel_size_real=18
   pixel_size=0.383   # in microns
   radii = np.arange(1,800)

#positions = [(512, 512)]
name_column='aperture_sum_'
#color=['red','blue','green','magenta','orange','black']

cmap = plt.get_cmap('gist_rainbow')
colors = ['grey','black']
colors.extend([cmap(i) for i in np.linspace(0, 1, len(psf_file[2:-2])-1)])

total_aperture_sums=[]

for psf in psf_file[1:-2]:

   data=fits.getdata('%s.fits' % psf)

   print (np.unravel_index(data.argmax(), data.shape))
   center=np.unravel_index(data.argmax(), data.shape)

   apertures = [CircularAperture(center, r=r) for r in radii]
   phot_table = aperture_photometry(data, apertures,method='exact')#,method='subpixel', subpixels=5)

   radius_inf=radii[-1]
   aperture_inf = CircularAperture(center, r=radius_inf)
   phot_table_inf = aperture_photometry(data, aperture_inf,method='exact')#,method='subpixel',subpixels=5)
   print (phot_table_inf)
   aperture_sums=[]
   for i in range(len(radii)):
       aperture_sums.append(phot_table[name_column+str(i)]/phot_table_inf['aperture_sum'])

   total_aperture_sums.append(aperture_sums)

for i,aper_sum in enumerate(total_aperture_sums):
   #plt.plot(radii*pixel_size,aper_sum,label=psf_file[-1].split()[i],color=color[i])
   plt.plot(radii*pixel_size/pixel_size_real,aper_sum,label=psf_file[-2].split()[i],color=colors[i])
#plt.axvline(x=7.5, linewidth=1, linestyle='--',color='red')
#plt.axvline(x=15, linewidth=1, linestyle='--',color='red')
#plt.axvline(x=22.5, linewidth=1, linestyle='--',color='red')
#plt.xlabel(r'Radius $\mu$m')
plt.xlabel(r'Radius around centroid (in Pixels)',fontsize=14)
plt.ylabel('Enclosed Energy',fontsize=14)
plt.xlim(0,10)
plt.ylim(0,1)
plt.xticks(np.arange(0,10.9,1))
plt.yticks(np.arange(0,1.05,0.1))
plt.legend(loc='lower right',prop={'size':10})
plt.grid(True)
plt.title('%s' % psf_file[-1])
plt.savefig('ee.png')
plt.show()

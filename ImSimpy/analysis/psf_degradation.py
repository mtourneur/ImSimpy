import os
import sys
import subprocess as sp
import numpy as np
from astropy.io import ascii
from astropy.table import Table,Column
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def  extract_degradation(input_file,percentage):
   """  Extractr the degradation for the desired completeness percentage. Percentage between 0 and 1  """

   completeness_cat=ascii.read('%s.txt' % input_file)

   nb_of_sources=sum(completeness_cat['col2'])
   nb_of_sources_det=sum(completeness_cat['col3'])
   #cdf_sources=np.zeros(len(completeness_cat['col2']))
   cdf_det=np.zeros(len(completeness_cat['col3']))
   for i in range(len(completeness_cat['col3'])):
       #cdf_sources[i]+=completeness_cat['col2'][i]
       if i == 0: cdf_det[i]=completeness_cat['col3'][i]
       cdf_det[i]=cdf_det[i-1]+completeness_cat['col3'][i]
   cdf_det/=nb_of_sources_det
   f = interp1d(cdf_det,completeness_cat['col1'],kind='linear')
   
   return f(percentage)


def find_mag(input_dir,bands,locs,seeings,zoom,percentage):

   for band in bands:
       input_cat_name='../data/catalog/input_sources_'+band
       for loc in locs:
           for seeing in seeings:
               result_file=input_dir+'completeness_%s_%s_s%s_p%s'  % (band,loc,seeing,str(int(round((percentage)*100))))
               res=[]
  
               if zoom[0] != None:
                   for z in zoom:
                       completeness_file=input_dir+'completeness_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                           
                       # extract the degradation for the desired completeness percentage
                       res.append([int(round((z-1.0)*100)),extract_degradation(completeness_file,percentage)])
               res=np.array(res)

               ascii.write(Table([res.T[0],res.T[1]],names=('degradation','mag')),'%s.txt' % result_file)


def plot_degradation_mag(input_dir,bands,locs,seeings,percentage):
   color =['blue','green','red'] # for the 3 seeings
   for band in bands:
       for loc in locs:
           plot_name='results/plots/mag_degradation_%s_%s_p%s'  % (band,loc,str(int(round((percentage)*100))))
           plt.clf()
           i=0
           for seeing in seeings:

               result_file=input_dir+'completeness_%s_%s_s%s_p%s'  % (band,loc,seeing,str(int(round((percentage)*100))))
               res=ascii.read('%s.txt' % result_file)
               plt.plot(res['degradation'],res['mag'],marker='o',color=color[i],label='seeing %s' % str(float(seeing)/10),lw=0.5)
               if seeing == seeings[0]: magmax = max(res['mag'])
               if seeing == seeings[-1]: magmin = min(res['mag'])
               #plt.plot(res['degradation'],res['mag'],lw='0.5',color=color[i])
               i+=1
           plt.xlabel('PSF degradation (in %)')
           plt.ylabel('Magnitude (AB)')
           plt.title('Magnitude reached for %s%% completeness / %s %s ' % (str(int(round((percentage)*100))),band,loc))
           plt.minorticks_on()
           plt.xlim(0,310)
           plt.ylim(magmin-0.3,magmax+0.3)
           plt.legend()
           plt.grid(b=True, which='major', linestyle='-')
           plt.grid(b=True, which='minor', linestyle='-',lw=0.3,color='grey')
           plt.savefig('%s.png' % plot_name)

def plot_completeness(input_dir,bands,locs,seeings,zoom):
   cmap = plt.get_cmap('gist_rainbow')
   colors=[cmap(i) for i in np.linspace(0, 1, len(zoom))]
   for band in bands:
       for loc in locs:
           for seeing in seeings:
               plot_name='results/plots/completeness_%s_%s_s%s_p%s'  % (band,loc,seeing,str(int(round((percentage)*100))))
               plt.clf()
               if zoom[0] != None:
                   i=0
                   for z in zoom:
                       completeness_file=input_dir+'completeness_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))

                       completeness_cat=ascii.read('%s.txt' % completeness_file)
                       
                       eff=completeness_cat['col3']/completeness_cat['col2']

                       mag=completeness_cat['col1']
                       plt.plot(mag,eff,marker='o',color=colors[i],label='%s%%' % str(int(round((z-1.0)*100))))
                       i+=1
                   plt.xlabel('Magnitude (AB)')
                   plt.ylabel('Detection efficiency')
                   plt.title('Detection efficiency / %s %s /seeing %s" ' % (band,loc,float(seeing)/10))
                   plt.minorticks_on()
                   #plt.xlim(21,22.5)
                   plt.ylim(0,1)
                   plt.legend(loc='lower right')
                   plt.grid(b=True, which='major', linestyle='-')
                   plt.grid(b=True, which='minor', linestyle='-',lw=0.3,color='grey')
                   plt.savefig('%s.png' % plot_name)

def max_degradation(input_dir,bands,locs,seeings,percentage,delta_mag):
   """ returns the maximum degradation for a delta mag """
   for band in bands:
       fileres=[]
       filename='delta_mag_%s_p%s_delta%s.txt' % (band,str(int(round(percentage*100))),str(int(round(delta_mag*100))))
       for loc in locs:
           for seeing in seeings:
               result_file=input_dir+'completeness_%s_%s_s%s_p%s'  % (band,loc,seeing,str(int(round((percentage)*100))))
               res=ascii.read('%s.txt' % result_file)
               f=interp1d(res['mag'],res['degradation'],kind='linear')
               if res['mag'][0]-delta_mag < res['mag'][-1]: degrad='>300'
               else: degrad = str(f(res['mag'][0]-delta_mag))
               fileres.append([band,loc,str(float(seeing)/10),degrad])
       np.savetxt(filename,fileres,fmt="%s %s %s %s")
if __name__ == '__main__':

   # Parameters for ds9 regions
   input_dir='results/psf_degrad_30s/'
   output_dir='results/psf_degrad_300s/'

   bands=['g','z','H']
   locs=['center','edge']
   seeings=['04','07','10']            # FWHM in arcsecond  between 0.025 and 5.5"
   zoom=[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,3.0,4.0]  #[None]
   percentage=0.8   
   delta_mag=0.7

   plot_completeness(input_dir,bands,locs,seeings,zoom)
   find_mag(input_dir,bands,locs,seeings,zoom,percentage)
   plot_degradation_mag(input_dir,bands,locs,seeings,percentage)
   max_degradation(input_dir,bands,locs,seeings,percentage,delta_mag)


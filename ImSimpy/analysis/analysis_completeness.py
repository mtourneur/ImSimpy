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


def find_mag(bands,locs,seeings,zoom,percentage):

   for band in bands:
       input_cat_name='../data/catalog/input_sources_'+band
       for loc in locs:
           for seeing in seeings:
               result_file='results/completeness_%s_%s_s%s_p%s'  % (band,loc,seeing,str(int(round((percentage)*100))))
               res=[]
  
               if zoom[0] != None:
                   for z in zoom:
                       if z == 1:
                           completeness_file='results/completeness_%s_%s_s%s'  % (band,loc,seeing)
                       else:
                           completeness_file='results/completeness_%s_%s_s%s_z%s'  % (band,loc,seeing,str(int(round((z-1.0)*100))))
                           
                       # extract the degradation for the desired completeness percentage
                       res.append([int(round((z-1.0)*100)),extract_degradation(completeness_file,percentage)])
               res=np.array(res)

               ascii.write(Table([res.T[0],res.T[1]],names=('degradation','mag')),'%s.txt' % result_file)


def plot_degradation_mag(bands,locs,seeings,percentage):

   for band in bands:
       input_cat_name='../data/catalog/input_sources_'+band
       for loc in locs:
           plot_name='results/plots/mag_degradation_%s_%s_p%s'  % (band,loc,str(int(round((percentage)*100))))
           for seeing in seeings:

               result_file='results/completeness_%s_%s_s%s_p%s'  % (band,loc,seeing,str(int(round((percentage)*100))))
               res=ascii.read('%s.txt' % result_file)
               plt.plot(res['degradation'],res['mag'],label='seeing %s' % str(float(seeing)/10))       

       plt.savefig('%s.png' % plot_name)

if __name__ == '__main__':

   # Parameters for ds9 regions
   output_dir='results/'

   bands=['z']
   locs=['center']
   seeings=['04','07','10']            # FWHM in arcsecond  between 0.025 and 5.5"
   zoom=[1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,3.0,4.0]  #[None]
   percentage=0.8   

   find_mag(bands,locs,seeings,zoom,percentage)
   plot_degradation_mag(bands,locs,seeings,percentage)


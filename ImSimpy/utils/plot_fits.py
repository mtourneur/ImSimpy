from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pyregion

inst_only=5
seeing='07'
zoom=['0','10','20','30','40','50','60','70','80','90','100','200','300']
band='z'
loc='edge'

if band=='g':
   pixsize=0.187
   pixsizereal=15
   res='0187'
   width=np.ceil(pixsizereal/pixsize)
elif band == 'z':
   pixsize=0.338
   pixsizereal=15
   res='0338'
   width=np.ceil(pixsizereal/pixsize)
elif band == 'H':
   pixsize=0.383
   pixsizereal=18
   res='0383'
   width=np.ceil(pixsizereal/pixsize)




"""
fig = plt.figure(1)
ax = plt.subplot(111)
plt.imshow(data,origin='lower',cmap='jet',vmin=0,vmax=1)
plt.colorbar()
patch_list, text_list = r.get_mpl_patches_texts()
for p in patch_list:
    ax.add_patch(p)
for t in text_list:
    ax.add_artist(t)

plt.axis('off')

#plt.tick_labels.hide()
#plt.axis_labels.hide()

plt.savefig('test2')
plt.show()
"""

if inst_only==1:

   fig, axs = plt.subplots(3,2, figsize=(13,6), sharex=True)
   pars = ['g center', 'g edge', 'z center', 'z edge', 'H center', 'H edge']
   print (pars)
   for i,p in enumerate(pars):
       print (i)
       if i ==0:
           data=fits.getdata('zemax/PSF_g_center.fits')
           center=[513,513]
           width=np.ceil(15/0.187)
       elif i ==1:
           data=fits.getdata('zemax/PSF_g_edge.fits')       
           center=[513,513]
           width=np.ceil(15/0.187)
       elif i ==2:
           data=fits.getdata('zemax/PSF_z_center.fits')
           center=[513,513]
           width=np.ceil(15/0.338)
       elif i==3:
           data=fits.getdata('zemax/PSF_z_edge.fits')
           center=[513,513]
           width=np.ceil(15/0.338)
       elif i==4:
           data=fits.getdata('zemax/PSF_H_center.fits')
           center=[513,513]
           width=np.ceil(18/0.383)
       elif i==5:
           data=fits.getdata('zemax/PSF_H_edge.fits')
           center=[513,513]
           width=np.ceil(18/0.383)
       data/=np.max(data)
       center=np.unravel_index(data[300:700,300:700].argmax(), data[300:700,300:700].shape)
       #center=[513,513]

       region_string = """
       # Region file format: DS9 version 4.1
       global color=black dashlist=8 3 width=1
       image 
       box(%d,%d,%d,%d,0)
       """ % (center[1],center[0],width,width)

       r = pyregion.parse(region_string)
       m = r.get_mask(shape=data.shape)
       patch_list, text_list = r.get_mpl_patches_texts()


       im=axs.flat[i].imshow(data[300:700,300:700],origin='lower',cmap='jet',vmin=0,vmax=1)
       axs.flat[i].text(0.95,0.85, p, transform=axs.flat[i].transAxes, ha='right',color='white',fontsize=11)
       axs.flat[i].axis('off')
       for p in patch_list:
           axs.flat[i].add_patch(p)
       for t in text_list:
           axs.flat[i].add_artist(t)

   #plt.suptitle('%s band %s / seeing: %.1f' % (band,loc,float(seeing)/10), fontsize=14)
   # Make an axis for the colorbar on the right side
   cax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
   fig.colorbar(im,cax=cax)
   #fig.tight_layout()
   plt.savefig('multiplot_inst')
   plt.show()

elif inst_only==2:

   fig, axs = plt.subplots(3,2, figsize=(13,6), sharex=True)
   pars = ['g center', 'g edge', 'z center', 'z edge', 'H center', 'H edge']
   print (pars)
   for i,p in enumerate(pars):
       print (i)
       if i ==0:
           data=fits.getdata('oversampled_psf/zemax_moffat_1024_g_center_s07.fits')
           x1,y1=100,900
           print (x1,y1)
           center=[513,513]
           width=np.ceil(15/0.187)
       elif i ==1:
           data=fits.getdata('oversampled_psf/zemax_moffat_1024_g_edge_s07.fits')
           x1,y1=100,900
           center=[513,513]
           width=np.ceil(15/0.187)
       elif i ==2:
           data=fits.getdata('oversampled_psf/zemax_moffat_1024_z_center_s07.fits')
           x1,y1=100,900
           center=[513,513]
           width=np.ceil(15/0.338)
       elif i==3:
           data=fits.getdata('oversampled_psf/zemax_moffat_1024_z_edge_s07.fits')
           x1,y1=100,900
           center=[513,513]
           width=np.ceil(15/0.338)
       elif i==4:
           data=fits.getdata('oversampled_psf/zemax_moffat_1024_H_center_s07.fits')
           x1,y1=100,900
           center=[513,513]
           width=np.ceil(18/0.383)
       elif i==5:
           data=fits.getdata('oversampled_psf/zemax_moffat_1024_H_edge_s07.fits')
           x1,y1=100,900
           center=[513,513]
           width=np.ceil(18/0.383)
       data/=np.max(data)
       center=np.unravel_index(data[x1:y1,x1:y1].argmax(), data[x1:y1,x1:y1].shape)
       #center=[513,513]

       region_string = """
       # Region file format: DS9 version 4.1
       global color=black dashlist=8 3 width=1
       image 
       box(%d,%d,%d,%d,0)
       """ % (center[1],center[0],width,width)

       r = pyregion.parse(region_string)
       m = r.get_mask(shape=data.shape)
       patch_list, text_list = r.get_mpl_patches_texts()

       
       im=axs.flat[i].imshow(data[x1:y1,x1:y1],origin='lower',cmap='jet',vmin=0,vmax=1)
       axs.flat[i].text(0.95,0.85, p, transform=axs.flat[i].transAxes, ha='right',color='white',fontsize=11)
       axs.flat[i].axis('off')
       for p in patch_list:
           axs.flat[i].add_patch(p)
       for t in text_list:
           axs.flat[i].add_artist(t)

   #plt.suptitle('%s band %s / seeing: %.1f' % (band,loc,float(seeing)/10), fontsize=14)
   # Make an axis for the colorbar on the right side
   cax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
   fig.colorbar(im,cax=cax)
   #fig.tight_layout()
   plt.savefig('multiplot_inst_seeing07')
   plt.show()



elif inst_only ==3:
   fig, axs = plt.subplots(4,4, figsize=(13,6), sharex=True)
   pars = 'Zemax 10% 20% 30% 40% 50% 60% 70% 80% 90% 100% 200% 300%'.split()
   for i,p in enumerate(pars):
      print (i)
      data=fits.getdata('psf_degradation/zemax_1024_%s_%s_z%s.fits' % (band,loc,zoom[i]))
      data/=np.max(data)
      center=np.unravel_index(data.argmax(), data.shape)
      #center=[513,513]

      region_string = """
      # Region file format: DS9 version 4.1
      global color=black dashlist=8 3 width=1
      image 
      box(%d,%d,%d,%d,0)
      """ % (center[1],center[0],width,width)

      r = pyregion.parse(region_string)
      m = r.get_mask(shape=data.shape)
      patch_list, text_list = r.get_mpl_patches_texts()


      im=axs.flat[i].imshow(data,origin='lower',cmap='jet',vmin=0,vmax=1)
      axs.flat[i].text(0.95,0.85, p, transform=axs.flat[i].transAxes, ha='right',color='white',fontsize=11)
      axs.flat[i].axis('off')
      for p in patch_list:
          axs.flat[i].add_patch(p)
      for t in text_list:
          axs.flat[i].add_artist(t)

   fig.delaxes(axs.flat[13])
   fig.delaxes(axs.flat[14])
   fig.delaxes(axs.flat[15])


   plt.suptitle('%s band %s' % (band,loc), fontsize=14)
   # Make an axis for the colorbar on the right side
   cax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
   fig.colorbar(im,cax=cax)
   #fig.tight_layout()
   plt.savefig('multiplot_%s_%s' % (band,loc))
   plt.show()


else:

   fig, axs = plt.subplots(4,4, figsize=(13,6), sharex=True)
   pars = '0% 10% 20% 30% 40% 50% 60% 70% 80% 90% 100% 200% 300% Zemax Moffat'.split()
   for i,p in enumerate(pars):
      print (i)
      if i ==13:
          data=fits.getdata('psf_degradation/zemax_1024_%s_%s_z%s.fits' % (band,loc,0))
      elif i==14:
          data=fits.getdata('atmosphere/moffat_1024_%s_s%s.fits' % (band,seeing))
      else:
          data=fits.getdata('psf_degradation/zemax_moffat_1024_%s_%s_s%s_z%s.fits' % (band,loc,seeing,zoom[i]))
      data/=np.max(data)
      center=np.unravel_index(data.argmax(), data.shape)
      #center=[513,513]
   
      region_string = """
      # Region file format: DS9 version 4.1
      global color=black dashlist=8 3 width=1
      image 
      box(%d,%d,%d,%d,0)
      """ % (center[1],center[0],width,width)
   
      r = pyregion.parse(region_string)
      m = r.get_mask(shape=data.shape)
      patch_list, text_list = r.get_mpl_patches_texts()
   
   
      im=axs.flat[i].imshow(data,origin='lower',cmap='jet',vmin=0,vmax=1)
      axs.flat[i].text(0.95,0.85, p, transform=axs.flat[i].transAxes, ha='right',color='white',fontsize=11)
      axs.flat[i].axis('off')
      for p in patch_list:
          axs.flat[i].add_patch(p)
      for t in text_list:
          axs.flat[i].add_artist(t)
   
   fig.delaxes(axs.flat[15])
   plt.suptitle('%s band %s / seeing: %.1f' % (band,loc,float(seeing)/10), fontsize=14)
   # Make an axis for the colorbar on the right side
   cax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
   fig.colorbar(im,cax=cax)
   #fig.tight_layout()
   plt.savefig('multiplot_%s_%s_s%s' % (band,loc,seeing))
   plt.show()


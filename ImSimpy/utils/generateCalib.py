"""
Generates Calibration Images 

"""
import os,datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
from astropy.io import fits

__version__ = 1.0


def Vignetting(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Vignetting/Vignetting_vis.fits',Type='quad', floor=1, xsize=4096, ysize=4096,create3Dplot=False,rstride=100,cstride=100):
    """

    """
    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if Type == 'quad':
        # generate random data
        x = (np.random.random(xsize) - 0.5)
        y = (np.random.random(ysize) - 0.5)
        xx, yy = np.meshgrid(np.linspace(x.min(), x.max(), xsize), np.linspace(y.min(), y.max(), ysize))
        surface = (-1.4*xx*xx - 1.6*yy*yy - 1.5*xx*yy)*floor*0.09 + floor        #about 9-10 per cent range
        #cutout extra
        surface = surface[:ysize, :xsize]

        if create3Dplot:
            x, y = np.meshgrid(np.arange(0, xsize, 1), np.arange(0, ysize, 1))

            #plot 3D
            fig = plt.figure()
            plt.title('VIS Flat Fielding: Idealised Calibration Unit Flux')
            ax = Axes3D(fig)
            ax.plot_surface(x, y, surface, rstride=rstride, cstride=cstride, alpha=0.5, cmap=cm.jet)
            ax.set_zlabel('electrons')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            plt.savefig('%s.pdf' % 'Vignetting')
            plt.close()

    #save to file
    hdu=fits.PrimaryHDU()

    #Normalised surface
    #surface/=np.max(surface)
    hdu.data=surface.astype(np.float32)

    if Type == 'constant': hdu.header.set('VIGNETT', 'Vignetting image  with quadratic formula')

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)


def Offset(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Offset/Offset_vis.fits',Type='constant', floor=1000, xsize=4096, ysize=4096):
    """

    """
    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if Type == 'constant':
        offset = np.ones([xsize,ysize])*floor

    """
    #plot 3D
    fig = plt.figure()
    plt.title('Bias offset image')
    ax = Axes3D(fig)
    ax.plot_surface(x, y, bias, rstride=100, cstride=100, alpha=0.5, cmap=cm.jet)
    ax.set_zlabel('ADU')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.savefig('bias.pdf')
    plt.close()
    """
    #save to file
    hdu=fits.PrimaryHDU()
    hdu.data=offset.astype(np.float32)

    if Type == 'constant': hdu.header.set('OFFSET', 'OFFSET image with constant offset of %d ADU' % int(floor))

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)

def GainMap(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/GainMap/Gain_vis.fits',Type='constant', mean=1, std=0.05,Nampl=32,xsize=4096, ysize=4096):
    """

    """
    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if Type == 'constant':
        gain_map = np.ones([xsize,ysize])*mean
    elif Type == 'random':
        Nampl=int(Nampl)
        gain = np.random.normal(loc=mean,scale=std,size=Nampl)
        gain_map=np.zeros([xsize,ysize])
        for i in range(Nampl):
            gain_map[:,int(xsize/Nampl)*i:int(xsize/Nampl)*(i+1)]=gain[i]

    """
    #plot 3D
    fig = plt.figure()
    plt.title('Bias offset image')
    ax = Axes3D(fig)
    ax.plot_surface(x, y, bias, rstride=100, cstride=100, alpha=0.5, cmap=cm.jet)
    ax.set_zlabel('ADU')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.savefig('bias.pdf')
    plt.close()
    """
    #save to file
    hdu=fits.PrimaryHDU()
    hdu.data=gain_map.astype(np.float32)

    if Type == 'constant': hdu.header.set('GAINMAP', 'Gain image with constant gain of %f' % mean)

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)


def DarkCurrent(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/DarkCurrent/DarkCurrent_vis.fits',Type='gaussian', mean=0.006,xsize=4096, ysize=4096,texp=20):
    """

    """
    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if Type == 'gaussian':
        DC = np.random.normal(loc=mean*texp,scale=np.sqrt(mean*texp),size=[xsize,ysize])

    #set negative values to zero
    DC[DC<0]=0
    #save to file
    hdu=fits.PrimaryHDU()
    hdu.data=DC.astype(np.float32)

    if Type == 'constant': hdu.header.set('DCMAP', 'Dark current with gaussian mean of %f' % mean)

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)

def Bias(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Bias/Bias_vis.fits',BiasType='constant',numdata=4096, floor=1000, xsize=4096, ysize=4096):
    """

    """
    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if BiasType == 'constant':
        bias = np.ones([xsize,ysize])*floor

    """
    #plot 3D
    fig = plt.figure()
    plt.title('Bias offset image')
    ax = Axes3D(fig)
    ax.plot_surface(x, y, bias, rstride=100, cstride=100, alpha=0.5, cmap=cm.jet)
    ax.set_zlabel('ADU')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.savefig('bias.pdf')
    plt.close()
    """
    #save to file
    hdu=fits.PrimaryHDU()
    hdu.data=bias.astype(np.float32)

    if BiasType == 'constant': hdu.header.set('BIAS', 'BIAS with constant offset of %d ADU' % int(floor))

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)


def Flat(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/FlatField/FlatField_vis.fits',FlatFieldType='quad',numdata=4096, floor=1e5, xsize=4096, ysize=4096,create3Dplot=False,rstride=100,cstride=100):
    """
      an idealised flat field surface representing the calibration unit flux output
    """
    # Create directory if not existing
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    if FlatFieldType == 'quad':
        # generate random data
        x = (np.random.random(numdata) - 0.5)
        y = (np.random.random(numdata) - 0.5)
        xx, yy = np.meshgrid(np.linspace(x.min(), x.max(), xsize), np.linspace(y.min(), y.max(), ysize))
        surface = (-1.4*xx*xx - 1.6*yy*yy - 1.5*xx*yy)*floor*0.09 + floor        #about 9-10 per cent range
        #cutout extra
        surface = surface[:ysize, :xsize]

        if create3Dplot:
            x, y = np.meshgrid(np.arange(0, xsize, 1), np.arange(0, ysize, 1))

            #plot 3D
            fig = plt.figure()
            plt.title('VIS Flat Fielding: Idealised Calibration Unit Flux')
            ax = Axes3D(fig)
            ax.plot_surface(x, y, surface, rstride=rstride, cstride=cstride, alpha=0.5, cmap=cm.jet)
            ax.set_zlabel('electrons')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            plt.savefig('%.pdf' % 'Flatfield')
            plt.close()

    #save to file
    hdu=fits.PrimaryHDU()

    #Normalised surface
    #surface/=np.max(surface)
    hdu.data=surface.astype(np.float32)

    if FlatFieldType == 'constant': hdu.header.set('FLATFIEL', 'Flatfield with quadratic formula')

    hdu.header.set('DATE-OBS', datetime.datetime.isoformat(datetime.datetime.now()))
    hdu.header.set('INSTRUME', 'ImSimpy_%s' % str(__version__))

    hdu.header.add_history('If questions, please contact David Corre (david.corre at lam.fr).')
    hdu.header.add_history('Created by ImSimpy (version=%.2f) at %s' % (__version__,
                                                                           datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.verify('fix')

    hdu.writeto(filename,overwrite=True)

    
def makeCalibs(output_dir='Calibration/COLIBRI/',suffix=['_vis','_nir'],xsize=[4096,2048],ysize=[4096,2048],offset=[1000,500],Type_vig=['quad','quad'],Type_gain=['constant','random'],mean_gain=[1.5,1.3],std_gain=[None,0.05],Nampl=[1,32]):
    # Make Offset, Vignetting and Gain map
    from ImSimpy.utils.generateCalib import Offset, Vignetting, GainMap
    
    for i in range(len(suffix)):
        Offset(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Offset/%s/Offset%s.fits' % (output_dir,suffix[i]),xsize=xsize[i], ysize=ysize[i],floor=offset[i])

        Vignetting(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Vignetting/%s/Vignetting%s.fits' % (output_dir,suffix[i]),Type=Type_vig[i],floor=1,xsize=xsize[i], ysize=ysize[i])

        GainMap(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/GainMap/%s/Gain%s.fits' % (output_dir,suffix[i]),Type=Type_gain[i],mean=mean_gain[i],std=std_gain[i],Nampl=Nampl[i],xsize=xsize[i],ysize=ysize[i])
    

    
def makeFlats(configFile=os.getenv('ImSimpy_DIR')+'/ImSimpy/configFiles/Flat.hjson',name_telescope='colibri',output_dir='Calibration/COLIBRI/',bands=['g','r','i','z','J','H'],Texp=30,Nexp=10,suffix=['_vis','_nir']):

    from ImSimpy.ImSimpy import ImageSimulator
    
    
    GFT_IS=ImageSimulator(configFile=configFile,name_telescope=name_telescope)
    #Read the configfile
    GFT_IS.readConfigs()

    GFT_IS.config['verbose']=False
    GFT_IS.config['etc_plot']=False
    GFT_IS.config['disp']=False

    for N in range(Nexp):
        N+=1    # Python array indexing starting at 0
        for band in bands:
   
            if band in ['gri','g','r','i','blue','B','V','R','I']:
                GFT_IS.config['channel']='VIS1'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[0])


            elif band in ['zy','z','y','red']:
                GFT_IS.config['channel']='VIS2'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[0])


            elif band in ['J','H']:
                GFT_IS.config['channel']='NIR'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[1])

            

            GFT_IS.config['exptime']=Texp
            GFT_IS.config['filter_band']=band
            GFT_IS.config['output']='%s/Flat_%s_%d.fits' % (output_dir,band,N)
            GFT_IS.simulate('data')
            
def makeDarks(configFile=os.getenv('ImSimpy_DIR')+'/ImSimpy/configFiles/Dark.hjson',name_telescope='colibri',output_dir='Calibration/COLIBRI/',bands=['g','r','i','z','J','H'],Texp=30,Nexp=10,suffix=['_vis','_nir']):

    from ImSimpy.ImSimpy import ImageSimulator
    
    
    GFT_IS=ImageSimulator(configFile=configFile,name_telescope=name_telescope)
    #Read the configfile
    GFT_IS.readConfigs()

    GFT_IS.config['verbose']=False
    GFT_IS.config['etc_plot']=False
    GFT_IS.config['disp']=False

    for N in range(Nexp):
        N+=1    # Python array indexing starting at 0
        for band in bands:
   
            if band in ['gri','g','r','i','blue','B','V','R','I']:
                GFT_IS.config['channel']='VIS1'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[0])


            elif band in ['zy','z','y','red']:
                GFT_IS.config['channel']='VIS2'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[0])


            elif band in ['J','H']:
                GFT_IS.config['channel']='NIR'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[1])

            

            GFT_IS.config['exptime']=Texp
            GFT_IS.config['filter_band']=band
            GFT_IS.config['output']='%s/Dark_%s_%ss_%d.fits' % (output_dir,band,Texp,N)
            GFT_IS.simulate('data')
            
def makeBias(configFile=os.getenv('ImSimpy_DIR')+'/ImSimpy/configFiles/Bias.hjson',name_telescope='colibri',output_dir='Calibration/COLIBRI/',bands=['g','r','i','z','J','H'],Texp=0.1,Nexp=10,suffix=['_vis','_nir']):

    from ImSimpy.ImSimpy import ImageSimulator
    
    
    GFT_IS=ImageSimulator(configFile=configFile,name_telescope=name_telescope)
    #Read the configfile
    GFT_IS.readConfigs()

    GFT_IS.config['verbose']=False
    GFT_IS.config['etc_plot']=False
    GFT_IS.config['disp']=False

    for N in range(Nexp):
        N+=1    # Python array indexing starting at 0
        for band in bands:
   
            if band in ['gri','g','r','i','blue','B','V','R','I']:
                GFT_IS.config['channel']='VIS1'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[0])


            elif band in ['zy','z','y','red']:
                GFT_IS.config['channel']='VIS2'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[0])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[0])


            elif band in ['J','H']:
                GFT_IS.config['channel']='NIR'
                GFT_IS.config['GainMapFile']='%s/Gain%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['VignettingFile']='%s/Vignetting%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['OffsetFile']='%s/Offset%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['DeadPixFile']='%s/DeadPixs%s.fits' % (output_dir,suffix[1])
                GFT_IS.config['HotPixFile']='%s/HotPixs%s.fits' % (output_dir,suffix[1])

            

            GFT_IS.config['exptime']=Texp
            GFT_IS.config['filter_band']=band
            GFT_IS.config['output']='%s/Bias_%s_%d.fits' % (output_dir,band,N)
            GFT_IS.simulate('data')
    
if __name__ == '__main__':
    # Visible
    Vignetting(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Vignetting/Vignetting_vis.fits',floor=1,xsize=4096, ysize=4096)
    Offset(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Offset/Offset_vis.fits',xsize=4096, ysize=4096)
    GainMap(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/GainMap/Gain_vis.fits',Type='constant',mean=1.5,xsize=4096,ysize=4096)
    #DarkCurrent(filename=os.getenv('IS_DIR')+'/IS/data/DarkCurrent/DarkCurrent_vis.fits',Type='gaussian',mean=0.006,xsize=4096,ysize=4096)
    Flat(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/FlatField/FlatField_vis.fits',floor=1,xsize=4096, ysize=4096)
    Bias(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Bias/Bias_vis.fits',xsize=4096, ysize=4096)



    # NIR
    Vignetting(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Vignetting/Vignetting_nir.fits',floor=1,xsize=2048, ysize=2048)
    Offset(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Offset/Offset_nir.fits',xsize=2048, ysize=2048)
    GainMap(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/GainMap/Gain_nir.fits',Type='random',mean=1.3,std=0.05,Nampl=32,xsize=2048,ysize=2048)
    #DarkCurrent(filename=os.getenv('IS_DIR')+'/IS/data/DarkCurrent/DarkCurrent_nir.fits',Type='gaussian',mean=0.002,xsize=2048,ysize=2048)
    Flat(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/FlatField/FlatField_nir.fits',floor=1,xsize=2048, ysize=2048)
    Bias(filename=os.getenv('ImSimpy_DIR')+'/ImSimpy/data/Bias/Bias_nir.fits',xsize=2048, ysize=2048)


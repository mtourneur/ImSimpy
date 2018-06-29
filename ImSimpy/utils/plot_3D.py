from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from astropy.io import fits



data=fits.getdata('zemax/PSF_z_edge.fits')
#data=fits.getdata('oversampled_psf/zemax_moffat_128_z_edge_s04.fits')
#data=fits.getdata('psf_degradation/zemax_moffat_1024_z_edge_s04_z300.fits')
#data=fits.getdata('airy_ideal_psfs/airy_1024_g.fits')

if np.sum(data)!=1: data/=np.sum(data)

fig = plt.figure()
ax = fig.gca(projection='3d')

center=np.unravel_index(data.argmax(), data.shape)
width=50
print (center)
# Make data.
X = np.arange(center[0]-width,center[0]+width,1)
Y = np.arange(center[1]-width,center[1]+width,1)
X, Y = np.meshgrid(X, Y)
Z = data[center[0]-width:center[0]+width,center[1]-width:center[1]+width]
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.jet,rstride=1, cstride=1,
                       linewidth=0, antialiased=False,shade=False)

# Customize the z axis.
#ax.set_zlim(0, 1.)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

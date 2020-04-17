import pylab as pl
import numpy as np

a = np.array([[-200, 400]])

pl.figure(figsize=(5, 1))

img = pl.imshow(a, cmap='gist_earth')
pl.gca().set_visible(False)
cax = pl.axes([0, 0, 0.03, 2])
cbar = pl.colorbar(cax=cax, extend='both', ticks=[-200, -100, 0, 100, 200, 300, 400], orientation='vertical')
cbar.ax.set_yticklabels([-200, -100, 0, 100, 200, 300, 400])
cax.yaxis.set_ticks_position('left')
#cbar.set_label('velocity gradient')
plt.savefig("cbar_elevation.png", bbox_inches='tight', dpi=200)





import pylab as pl
import numpy as np

a = np.array([[0, 30]])

pl.figure(figsize=(5, 1))

img = pl.imshow(a, cmap='jet')
pl.gca().set_visible(False)
cax = pl.axes([0, 0, 0.03, 2])
cbar = pl.colorbar(cax=cax, extend='both', ticks=[0, 10, 20, 30], orientation='vertical')
cbar.ax.set_yticklabels([0, 10, 20, 30])
#cbar.set_label('velocity gradient')
plt.savefig("cbar_ice_velo.png", bbox_inches='tight', dpi=200)





import pylab as pl
import numpy as np

a = np.array([[1250, 2500]])

pl.figure(figsize=(5, 1))

img = pl.imshow(a, cmap='RdYlBu')
pl.gca().set_visible(False)
cax = pl.axes([0, 0, 0.03, 2])
cbar = pl.colorbar(cax=cax, extend='both', ticks=[1250, 1500, 1750, 2000, 2250, 2500], orientation='vertical')
cbar.ax.set_yticklabels([1250, 1500, 1750, 2000, 2250, 2500])
cax.yaxis.set_ticks_position('left')
#cbar.set_label('velocity gradient')
plt.savefig("cbar_thickness.png", bbox_inches='tight', dpi=200)





import pylab as pl
import numpy as np

a = np.array([[0, 0.2]])

pl.figure(figsize=(5, 1))

img = pl.imshow(a, cmap='viridis')
pl.gca().set_visible(False)
cax = pl.axes([0, 0, 0.03, 2])
cbar = pl.colorbar(cax=cax, extend='both', ticks=[0, 0.1, 0.2], orientation='vertical')
cbar.ax.set_yticklabels(['0.0', '0.1', '0.2'])
#cbar.set_label('velocity gradient')
plt.savefig("cbar_velograd.png", bbox_inches='tight', dpi=200)
#!/usr/bin/env python

#
# An example of using the weirdoutils to access AutoZ results
#

import numpy as np
import matplotlib.pyplot as plt
import weirdoutils as w


# We want to loop over all so this gets us the list of objects
info = w.load_info()

# The layout
nx_plot, ny_plot = 7, 4
dx = 1
dy = 2
# figsize = plt.figaspect(float(dy * ny_plot) / float(dx * nx_plot))
fig, axes = plt.subplots(nx_plot, ny_plot, figsize=(8, 15))
plt.subplots_adjust(hspace=0.2, wspace=0.05)
plt.setp(axes, xticks=[], yticks=[])

# Loop over it all.
for i, name in enumerate(info['Name']):

    a = w.load_autoz_results(name, bin=3, toget=['ZMAP'])

    nx, ny = a['ZMAP'].shape

    xc = np.floor(nx/2.).astype(int)
    yc = np.floor(ny/2.).astype(int)


    if name == 'SDSS024639-080525':
        yc, xc = 63, 91
    elif name == 'SDSS112934+032304':
        yc, xc = 43, 55
    elif name == 'SDSS024639-080525':
        yc, xc = 63, 91

    objz = np.median(a['ZMAP'][xc-2:xc+2, yc-2:yc+2])

    # The velocity map relative to systematic
    vmap = 3e5*(a['ZMAP']-objz)

    i_y = int(i/nx_plot)
    i_x = i - i_y*nx_plot

    # For display it seems that 
    axes[i_x, i_y].imshow(vmap, vmin=-140, vmax=140, cmap='RdBu')

    axes[i_x, i_y].set_title(name, horizontalalignment='center',
                            color='purple', fontsize=8)

fig.savefig('velocity-fields.png')
#plt.show()


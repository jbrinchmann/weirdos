#!/usr/bin/env python

#
# An example of using the weirdoutils to access data from a platefit fit. 
#

import numpy as np
import matplotlib.pyplot as plt
import weirdoutils as w
import sys

try:
    what = sys.argv[1]
    min = sys.argv[2]
    max = sys.argv[3]
except:
    print("""
Usage: 
    show_one_platefit_quantity.py <NAME> min max
    """)

# We want to loop over all so this gets us the list of objects
info = w.load_info()

# The layout
nx_plot, ny_plot = 7, 4

#dx = 1
#dy = 2
#figsize = plt.figaspect(float(dy * ny_plot) / float(dx * nx_plot))

fig, axes = plt.subplots(nx_plot, ny_plot, figsize=(8, 15))
plt.subplots_adjust(hspace=0.2, wspace=0.05)
plt.setp(axes, xticks=[], yticks=[])

# Loop over it all.
for i, name in enumerate(info['Name']):

    # Get the platefit result.
    d = w.load_platefit_results(name, what, bin=3, full_size=True)
    x = d[what]

    # Figure out the shape and the centre
    nx, ny = x.shape

    xc = np.floor(nx/2.).astype(int)
    yc = np.floor(ny/2.).astype(int)


    if name == 'SDSS024639-080525':
        yc, xc = 63, 91
    elif name == 'SDSS112934+032304':
        yc, xc = 43, 55
    elif name == 'SDSS024639-080525':
        yc, xc = 63, 91

    i_y = int(i/nx_plot)
    i_x = i - i_y*nx_plot

    # For display it seems that 
    im = axes[i_x, i_y].imshow(x, vmin=min, vmax=max)

    axes[i_x, i_y].set_title(name, horizontalalignment='center',
                            color='purple', fontsize=8)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
cbar_ax.set_title(what)

fig.savefig(what+'-distributions.png')
#plt.show()


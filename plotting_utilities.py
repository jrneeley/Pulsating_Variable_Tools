import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch, ImageNormalize
from astropy.io import fits
from matplotlib.patches import Circle

# Helper function to plot region of image
def plot_region(x, y, image, xall=[], yall=[], ext=0,
    axes=None, xoff=0, yoff=0, aperture=None, img_limits=[0,500]):

    image_data = fits.getdata(image, ext=ext)

    if axes is None:
        fig, ax = plt.subplots(1,1, figsize=(8,5))
    else:
        ax = axes

    norm1 = ImageNormalize(image_data, vmin=img_limits[0], vmax=img_limits[1],
        stretch=LogStretch())

    ax.imshow(image_data, cmap='gray', norm=norm1)

    #ax.set_aspect('equal')
    if aperture != None:
        ap = Circle((x-1, y-1), aperture, facecolor=None, edgecolor='red', fill=0)
        ax.add_patch(ap)
    if (len(xall) > 0) & (len(yall) > 0):
        x_all = xall - (xoff+1)
        y_all = yall - (yoff+1)
        ax.scatter(x_all, y_all, marker='x', color='green')
    ax.set_xlim(x-20, x+20)
    ax.set_ylim(y-20, y+20)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    if axes == None:
        return fig, ax
    #return fig, ax

# Helper function to plot LMC PL relations 
def plot_lmc_pl(axes=None, offset=0, period_cutoff=0, colors=True):

    if axes == None:
        fig1, ax1 = plt.subplots(1,1)

    else:
        ax1 = axes

    if colors == True:
        #colors = ['xkcd:sage', 'xkcd:gray', 'xkcd:pale purple', 'xkcd:rose',
        #    'xkcd:steel blue', 'xkcd:puce', 'xkcd:eggplant']
        colors = ['xkcd:magenta', 'xkcd:magenta', 'xkcd:aquamarine', 'xkcd:aquamarine',
            'xkcd:french blue', 'xkcd:goldenrod', 'xkcd:goldenrod']
    else:
        colors = ['xkcd:gray' for i in range(7)]

    al = 0.15
    # classical cepheid lines
    x_fo = np.array([-0.6, 0.8])
    x_fu = np.array([0.0, 2.1])
    if period_cutoff != 0:
        x_fo[1] = np.min([x_fo[1], np.log10(period_cutoff)])
        x_fu[1] = np.min([x_fu[1], np.log10(period_cutoff)])
    y_fo = -3.311*(x_fo-1.0) + 12.897 - 18.477 + offset
    y_fu = -2.912*(x_fu-1.0) + 13.741 - 18.477 + offset
    ax1.fill_between(x_fo, y_fo-0.16, y_fo+0.16, color=colors[0], alpha=al)
    ax1.fill_between(x_fu, y_fu-0.15, y_fu+0.15, color=colors[1], alpha=al)
    # anomalous cepheid lines
    x_fo = np.array([-0.4, 0.07])
    x_fu = np.array([-0.2, 0.37])
    if period_cutoff != 0:
        x_fo[1] = np.min([x_fo[1], np.log10(period_cutoff)])
        x_fu[1] = np.min([x_fu[1], np.log10(period_cutoff)])
    y_fo = -3.302*x_fo + 16.656 - 18.477 + offset
    y_fu = -2.962*x_fu + 17.368 - 18.477 + offset
    ax1.fill_between(x_fo, y_fo-0.16, y_fo+0.16, color=colors[2], alpha=al)
    ax1.fill_between(x_fu, y_fu-0.23, y_fu+0.23, color=colors[3], alpha=al)
    # type 2 cepheid line
    x_fu = np.array([-0.09, 1.8])
    if period_cutoff != 0:
        x_fu[1] = np.min([x_fu[1], np.log10(period_cutoff)])
    y_fu = -2.033*x_fu + 18.015 - 18.477 + offset
    ax1.fill_between(x_fu, y_fu-0.4, y_fu+0.4, color=colors[4], alpha=al)
    # RRL lines
    x_fo = np.array([-0.7, -0.3])
    x_fu = np.array([-0.6, 0.0])
    if period_cutoff != 0:
        x_fo[1] = np.min([x_fo[1], np.log10(period_cutoff)])
        x_fu[1] = np.min([x_fu[1], np.log10(period_cutoff)])
    y_fo = -2.014*x_fo + 17.743 - 18.477 + offset
    y_fu = -1.889*x_fu + 18.164 - 18.477 + offset
    ax1.fill_between(x_fu, y_fu-0.15, y_fu+0.15, color=colors[5], alpha=al)
    ax1.fill_between(x_fo, y_fo-0.16, y_fo+0.16, color=colors[6], alpha=al)

    if axes == None:
        ax1.set_xlabel('$\log P$')
        ax1.set_ylabel('I mag')
        ax1.invert_yaxis()
        ax2.set(xlabel='P [days]', ylabel='I amp')
        plt.show()

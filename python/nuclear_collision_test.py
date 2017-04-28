import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pylab import *

def loadEnergyDensity(t, i, nx, ny, nz, plotDir):
    hbarc=0.197327

    outputDir = root + '/cpu_vah_output'
    x, y, n, eRaw = np.loadtxt(outputDir + '/e_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    e_ideal = np.reshape(eRaw, (nx, ny, nz))*hbarc
    cmin = np.min(np.min(np.min(e_ideal)))
    cmax = np.max(np.max(np.max(e_ideal)))
    e_ideal_2p1=squeeze(e_ideal)

    outputDir = root + '/cpu_vah_output_pimunu'
    x, y, n, eRaw = np.loadtxt(outputDir + '/e_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    e_shear = np.reshape(eRaw, (nx, ny, nz))*hbarc
    #cmin = np.min(cmin,np.min(np.min(np.min(e_shear))))
    #cmax = np.max(np.max(np.max(np.max(e_shear))))
    e_shear_2p1=squeeze(e_shear)

    outputDir = root + '/cpu_vah_output_pimunu_Pi'
    x, y, n, eRaw = np.loadtxt(outputDir + '/e_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    e_shearBulk = np.reshape(eRaw, (nx, ny, nz))*hbarc
    #cmin = np.min(cmin,np.min(np.min(np.min(e_shearBulk))))
    #cmax = np.max(np.max(np.max(np.max(e_shearBulk))))
    e_shearBulk_2p1=squeeze(e_shearBulk)

    cmap='jet'

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16

    extent = [-lim, lim, -lim, lim]
    labelSize = 16

    plt.figure(i)
    plt.imshow(e_ideal_2p1, vmin=cmin, vmax=cmax, cmap=cmap, origin='lower', extent=extent)
    plt.title(r'$\tau=$ ' + str(t) + ' fm/c, ideal', fontsize=labelSize)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$y\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.colorbar()
    savefig(plotDir + '/e_ideal' + '{:.3f}'.format(t) + '.pdf', bbox_inches='tight')

    plt.figure(i+1)
    plt.imshow(e_shear_2p1, vmin=cmin, vmax=cmax, cmap=cmap, origin='lower', extent=extent)
    plt.title(r'$\tau=$ ' + str(t) + ' fm/c, shear', fontsize=labelSize)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$y\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.colorbar()
    savefig(plotDir + '/e_shear' + '{:.3f}'.format(t) + '.pdf', bbox_inches='tight')

    plt.figure(i+2)
    plt.imshow(e_shearBulk_2p1, vmin=cmin, vmax=cmax, cmap=cmap, origin='lower', extent=extent)
    plt.title(r'$\tau=$ ' + str(t) + ' fm/c, shear/bulk', fontsize=labelSize)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$y\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.colorbar()
    savefig(plotDir + '/e_shearBulk' + '{:.3f}'.format(t) + '.pdf', bbox_inches='tight')

if __name__ == '__main__':
    root = '/media/bazow/Data'
    plotDir = 'tests/figs/qgp/vah'

    nx = 161
    ny = 161
    nz = 1
    dx = 0.1
    dy = 0.1
    dz = 0.05
    t = 0.5
    tlist=[0.5, 2, 6, 10]

    lim = (nx - 1) / 2 * dx

    nt = len(tlist)

    k=0
    for i in range(0, nt):
        t = tlist[i]
        loadEnergyDensity(t, k, nx, ny, nz, plotDir)
        k += 3

    plt.show()

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pylab import *

if __name__ == '__main__':
    plotDir = 'tests/figs/qgp/vah/updated'

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

    hbarc=0.197327

    cmap='jet'

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16

    extent = [-lim, lim, -lim, lim]
    labelSize = 16

    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/'
    type = 'cpu-vah/'
    tag='mcGlb_2d_pimunu_isotropicIC_QCDEOS_etaOver20p2_reg0p1-1_theta1p1'

    dataDir = root + type + tag

    t=2
    x, y, n, eRaw = np.loadtxt(dataDir + '/e_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    e_shearBulk = np.reshape(eRaw, (nx, ny, nz))*hbarc
    e=squeeze(e_shearBulk)

    plt.figure(0)
    plt.imshow(e, cmap=cmap, origin='lower', extent=extent)
    plt.title(r'$\tau=$ ' + str(t) + ' fm/c, ideal', fontsize=labelSize)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$y\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.colorbar()
    savefig(plotDir + '/e_ideal' + '{:.3f}'.format(t) + '.pdf', bbox_inches='tight')

    t=6
    x, y, n, eRaw = np.loadtxt(dataDir + '/e_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    e_shearBulk = np.reshape(eRaw, (nx, ny, nz))*hbarc
    e=squeeze(e_shearBulk)

    plt.figure(1)
    plt.imshow(e, cmap=cmap, origin='lower', extent=extent)
    plt.title(r'$\tau=$ ' + str(t) + ' fm/c, shear', fontsize=labelSize)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$y\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.colorbar()
    savefig(plotDir + '/e_shear' + '{:.3f}'.format(t) + '.pdf', bbox_inches='tight')

    t=9.85
    x, y, n, eRaw = np.loadtxt(dataDir + '/e_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    e_shearBulk = np.reshape(eRaw, (nx, ny, nz))*hbarc
    e=squeeze(e_shearBulk)

    plt.figure(2)
    plt.imshow(e, cmap=cmap, vmin=0, vmax=0.7, origin='lower', extent=extent)
    plt.title(r'$\tau=$ ' + str(10) + ' fm/c, shear/bulk', fontsize=labelSize)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$y\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.colorbar()
    savefig(plotDir + '/e_shearBulk' + '{:.3f}'.format(10.0) + '.pdf', bbox_inches='tight')


    plt.show()

#!/usr/bin/env python3
from pylab import *
import StringIO


def loadVar(dir, var, t):
    buf = StringIO.StringIO()
    buf.write(dir+'/'+var+'_%.3f.dat' % (t))
    x, y, n, res = np.loadtxt(buf.getvalue(), unpack=True)
    return res

if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/cpu-vah'
    dataDir = root + '/nonconformal_bjorken_test'
    outputDir = 'tests/output/bjorken/nonconformal'
    
    t0 = 0.5
    tf = 20.0
    dt = 0.05
    nt = np.floor((tf - t0) / dt) + 1
    t = np.zeros((nt, 1))
    e = np.zeros((nt, 1))
    pl = np.zeros((nt, 1))
    pi = np.zeros((nt, 1))

    k = 0
    while k < nt:
        ti = t0 + k * dt
        t[k] = ti
        e[k] = loadVar(dataDir, 'e', ti)
        pl[k] = loadVar(dataDir, 'pl', ti)
        pi[k] = loadVar(dataDir, 'Pi', ti)
        k += 1

#    plt.figure(1)
#    plt.semilogy(t, e / e[0], linestyle='--', label=r'$\mathcal{E}/\mathcal{E}_{0}$')
#    plt.semilogy(t, pl / pl[0], 'r--', label=r'$\pi/\pi_{0}$')
#    plt.tight_layout(pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])

#    plt.show()

    np.savetxt(outputDir+'/e.dat', np.c_[t,e], fmt="%.16f")
    np.savetxt(outputDir+'/pl.dat', np.c_[t,pl], fmt="%.16f")
    np.savetxt(outputDir+'/Pi.dat', np.c_[t,pi], fmt="%.16f")

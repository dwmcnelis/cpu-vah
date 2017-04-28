#!/usr/bin/env python3

from scipy import integrate

import numpy
from matplotlib.pylab import *

import equation_of_state as eos
import specific_bulk_viscosity as zetas
#from plot_setup import plt

def zeta_zz_tc(a,b,pL):
    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7 = a6*a
    a8 = a7*a
    a9 = a8*a
    a10 = a9*a
    a11 = a10*a
    zetaLL0 = (-2.169274157232591e-6 + 0.0015007723408727923*a + 0.234244347953902*a2 + 2.956646117194832*a3 - 1.4546308078344168*a4 - 11.821812720171593*a5 +
     15.232146387176067*a6 - 5.145978753104648*a7)/(0.07331160572576197 + 3.120097137276907*a + 12.778215835624462*a2 - 31.764095787266974*a3 - 2.7994204846452098*a4 + 37.043802032024445*a5 -
     20.984106677521524*a6 + 2.539148721668404*a7)
    zetaLL2 = (0.00007482123591721489 - 0.0517637490482641*a - 8.079427938560013*a2 - 101.97917688859472*a3 + 50.172161480445254*a4 + 407.75277017721646*a5 -
     525.3802821796353*a6 + 177.4927576883693*a7)/   (30.34346422237839 + 1291.401068102149*a + 5288.883823777771*a2 - 13147.09924882401*a3 - 1158.6918230835379*a4 + 15332.384128302825*a5 - 8685.29447768419*a6 +
     1050.9506440540636*a7)
    Ihatm24 = 	T**4*(zetaLL0+36*b*zetaLL2)
    #return Ihatm24-3*pL
    return T**4*zetaLL0

def g0m2040DividedBym2_fun(a,b,e,p):
    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7 = a6*a
    a8 = a7*a
    a9 = a8*a
    a10 = a9*a
    a11 = a10*a
    meqHat2 = 36*b
    z2Logz = -2*(1-2/9*(e-3*p)/(b*(e+p))+0.02536*meqHat2)
    g0m2040DividedBym2ConstantTerm = (0.0002022034619596147 - 0.3465776971894224*a - 3669.180338013278*a10 - 0.46483928190644463*a2 + 305.8182449301262*a3 - 3094.0393881108607*a4 +
     13530.068145958274*a5 - 33008.96608710712*a6 + 47814.11155495329*a7 - 40754.13233671886*a8 + 18877.099618956072*a9)/   (0.3231726678495768 + 6.510257837715722*a - 146.6399433638468*a10 - 66.42749605402405*a2 + 278.2072846508577*a3 - 659.8881404723605*a4 +
     1009.6856692105162*a5 - 1129.2454150937733*a6 + 1136.3819011323233*a7 - 1029.3638914685152*a8 + 600.4581631072344*a9)
    g0m2040DividedBym2Zm2Term = (-5.428668780471906e-6 + 0.07087180484499124*a + 14430.23620485441*a10 - 3.6056344041371626*a2 - 75.98232257237731*a3 + 314.3389560219402*a4 -
     4860.077726215435*a5 + 28903.651438438494*a6 - 74655.13722955724*a7 + 95355.62148623986*a8 - 59409.22212071206*a9)/   (0.06525900311126431 + 1.7911514139670908*a + 284.5874952982027*a10 - 1.0324279865217276*a2 + 81.69891535142081*a3 - 259.5916250456495*a4 -
     59.79985417459323*a5 + 2029.823700195329*a6 - 4622.384249181224*a7 + 4479.931023062334*a8 - 1935.0883554632262*a9)
    g0m2040DividedBym2LogZTerm = (1.5260366926367532e-7 - 0.00024895619154304205*a - 16.08135726710075*a10 + 2.635895081202281*a11 - 0.0027342593291142404*a2 + 0.20899592202076162*a3 -
     2.1908540237701564*a4 + 10.836996467643049*a5 - 31.053876569731955*a6 + 55.16397490084364*a7 - 61.60172864783581*a8 + 42.08494163892432*a9)/   (0.0000605632612588605 + 0.0012923594925415354*a - 0.18379144108413883*a10 + 0.038614228316574885*a11 - 0.016304038454467058*a2 + 0.07975759323101335*a3 -
     0.21819092742589902*a4 + 0.37685396820303946*a5 - 0.45029351728017364*a6 + 0.44970998105242493*a7 - 0.44680890883509927*a8 + 0.3691002339464134*a9)
#    return (g0m2040DividedBym2LogZTerm * z2Logz + g0m2040DividedBym2Zm2Term) / meqHat2 + g0m2040DividedBym2ConstantTerm
    return g0m2040DividedBym2ConstantTerm

def g0m2020_fun(a,b,e,p):
    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7 = a6*a
    a8 = a7*a
    a9 = a8*a
    a10 = a9*a
    a11 = a10*a
    meqHat2 = 36*b
    z2Logz = -2*(1-2/9*(e-3*p)/(b*(e+p))+0.02536*meqHat2)
    g0m2020ConstantTerm = (-0.000011732877574172216 + 0.0008995982417255349*a + 0.7361757170685882*a2 + 8.855369493730423*a3 - 76.00475510029524*a4 + 261.7208147369771*a5 -
        479.7229720410445*a6 + 487.6088140168412*a7 - 260.7802309619955*a8 + 57.58861537142237*a9)/ (
        -0.0001445869624012341 + 0.03556951252765837*a + 3.2344009306364274*a2 + 12.514996006390815*a3 - 118.07971915906703*a4 + 347.32729764112844*a5 -
        483.55235351080324*a6 + 310.7815908976544*a7 - 67.06626623263872*a8 - 5.19154598317511*a9)
    g0m2020Z2Term = (-0.02907781069921067 - 13.745970067614193*a + 13997.84367889844*a10 - 325.2706635847172*a2 + 2604.8287335227815*a3 - 4049.1372422838062*a4 -
        16667.365141623744*a5 + 83158.20010908133*a6 - 154605.3533428809*a7 + 147412.71289551136*a8 - 71512.66891762966*a9)/ (
        0.9697656616566325 + 117.19185209632286*a + 7621.0452633237755*a10 + 842.3723825380051*a2 - 5926.237511293323*a3 + 14247.062719381864*a4 -
        16353.128301113404*a5 + 21213.272784815457*a6 - 46076.07167476786*a7 + 59110.92666130391*a8 - 34796.81744513102*a9)
    g0m2020Z2LogZTerm = (-0.0005534368375549849 - 0.3080130050685187*a + 630.3245252759478*a10 - 11.136413628061584*a2 + 34.63107256159966*a3 + 367.6376955018988*a4 -
                         2567.818354775122*a5 + 7178.293271719852*a6 - 10678.12443252271*a7 + 8758.211709132715*a8 - 3711.7083581405914*a9)/(
        0.009248889228618452 + 1.4181948328169438*a + 95.39051233960058*a10 + 21.503144137961534*a2 - 62.705284754794086*a3 - 83.17143447631213*a4 +
        693.6443740085989*a5 - 1134.2586687505138*a6 + 513.1520937342799*a7 + 331.751797652403*a8 - 376.7067817672707*a9)
    #return g0m2020ConstantTerm + g0m2020Z2LogZTerm * z2Logz + g0m2020Z2Term * meqHat2
    return g0m2020ConstantTerm

def zeta_z_tc(a,b):
    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7 = a6*a
    a8 = a7*a
    a9 = a8*a
    a10 = a9*a
    a11 = a10*a
    zetaLR = (1.6789989237359663e-6 - 0.00064490984848377*a - 0.08372751988817427*a2 - 0.5320705463072806*a3 + 2.7041802331042777*a4 - 4.064569383338656*a5 +
     2.5684272810257744*a6 - 0.5915938736175316*a7)/(-0.00019686722961534964 + 0.10476915315451184*a + 7.104762231638746*a2 + 30.346380785241006*a3 - 137.39710807489823*a4 + 162.51277077058876*a5 -
     70.87384079846204*a6 + 8.191285063886802*a7)
    return 36*b*T**4*zetaLR

def equations(t, y):
    delta_pipi = 1.33333
    tau_pipi = 1.42857
    delta_PiPi = 0.666667
    lambda_piPi = 1.2
    etabar = 0.2

    # Assign some variables for convenience of notation
    e = y[0]
    pL = y[1]
    pi = y[2]

    a = pL/e

    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7 = a6*a
    a8 = a7*a
    a9 = a8*a
    a10 = a9*a
    a11 = a10*a

    T = eos.effectiveTemperature(e)
    taupiInv = T / 5 / etabar
    p = eos.equilibriumPressure(e)
    cs2 = eos.speedOfSoundSquared(e)
    b = 1.0 / 3.0 - cs2
    if b==0: b=1e-16
    b2 = b * b
    zetabar = zetas.bulkViscosityToEntropyDensity(T)
    tauPiInv = 15 * b2 * T / zetabar

    zeta_zz = zeta_zz_tc(a,b,pL)

    g0m2040DividedBym2 = g0m2040DividedBym2_fun(a,b,e,p)
    beta_zPi = 3*g0m2040DividedBym2

    zeta_z = zeta_z_tc(a,b)

    beta_PiPi = 1-g0m2020_fun(a,b,e,p)

    # Output from ODE function must be a COLUMN vector, with n rows
    n = len(y)  # 2: implies we have two ODEs
    f = np.zeros((n, 1))
    f[0] = -(e + pL) / t
    f[1] = -(pL - p) * taupiInv + zeta_zz / t - beta_zPi * pi / t
    f[2] = -pi * tauPiInv - zeta_z / t - beta_PiPi * pi / t
    return f

if __name__ == '__main__':
    r = integrate.ode(equations).set_integrator('vode', method='bdf')

    t0 = 0.25
    t_final = 20.0
    delta_t = 0.05
    num_steps = np.floor((t_final - t0) / delta_t) + 1

    # Set initial condition
    T0 = 3.05
    ed = eos.equilibriumEnergyDensity(T0)
    T = eos.effectiveTemperature(ed)
    pinn = -2 / (3 * t0 * t0 * t0) * 0.2 * (ed + eos.equilibriumPressure(ed)) / T
    pizz0 = -t0 * t0 * pinn
    pi0 = 0
    pL0 = eos.equilibriumPressure(e)
    r.set_initial_value([ed, pL0, pi0], t0)

    t = np.zeros((num_steps, 1))
    e= np.zeros((num_steps, 1))
    pL = np.zeros((num_steps, 1))
    pi = np.zeros((num_steps, 1))
    t[0] = t0
    e[0] = ed
    pL[0] = pizz0
    pi[0] = pi0

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + delta_t)

        # Store the results to plot later
        t[k] = r.t
        e[k] = r.y[0]
        pL[k] = r.y[1]
        pi[k] = r.y[2]
        k += 1

    dataDir = 'tests/output/bjorken/nonconformal'

    tEst, eEst = np.loadtxt(dataDir+'/e.dat', unpack=True)
    tEst, plEst = np.loadtxt(dataDir+'/pl.dat', unpack=True)
    tEst, piEst = np.loadtxt(dataDir+'/Pi.dat', unpack=True)

    plotDir = 'tests/figs/bjorken/nonconformal'

    hbarc = 0.197327

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    labelSize = 16

    plt.figure(1)
    #ax = plt.semilogy(t, e / e[0], color='dodgerblue', alpha=0.45, linewidth=5, linestyle='-')
    #plt.semilogy(t, pL / pL[0], color='red', linewidth=5, alpha=0.25, linestyle='-')
    #ax = plt.plot(t, e / e[0], color='dodgerblue', alpha=0.45, linewidth=5, linestyle='-')
    #plt.plot(t, pL / pL[0], color='red', linewidth=5, alpha=0.25, linestyle='-')
    plt.plot(tEst, eEst / eEst[0], linestyle='--', label=r'$\mathcal{E}/\mathcal{E}_{0}$')
    plt.plot(tEst, 2*numpy.divide(plEst,eEst-plEst), 'r--', label=r'$\pi/\pi_{0}$')
    plt.xlabel(r'$\tau\,[\mathrm{fm/c}]$', fontsize=labelSize)
    plt.ylabel(r'${}\,\mathrm{{}{}}^{}$', fontsize=labelSize)
    plt.tight_layout(pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])
    plt.legend(loc='best', frameon=False, fontsize=labelSize)
    savefig(plotDir+'/nonconformalBjorkenFig.pdf', pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])
    #
    plt.figure(2)
#    plt.plot(t, pi*hbarc, color='dodgerblue', linewidth=5, alpha=0.25, linestyle='-', label='semi-analytic')
    plt.plot(tEst, piEst*hbarc, color=u'#1f77b4', linestyle='--', label='GPU-VH')
    plt.ylabel(r'$\Pi\,[\mathrm{GeV/fm}^{3}]$', fontsize=labelSize)
    plt.xlabel(r'$\tau\,[\mathrm{fm/c}]$', fontsize=labelSize)
    plt.legend(loc='best', frameon=False, fontsize=labelSize)
    plt.tight_layout(pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])
    savefig(plotDir+'/nonconformalBjorkenFig_Pi.pdf', pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])

    plt.show()

"""
Python script to generate a parameter space plot of
close encounters between a disk and a flyby object.

Author: Guang-Xing Li (SWIFAR), Xing Lu (SHAO)
Last update: Apr 5 2022
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from scipy.special import kn

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
rcParams.update({'font.size': 10})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = 'on'
rcParams['ytick.right'] = 'on'
rcParams['axes.xmargin'] = 0
rcParams['xtick.major.size'] = 4
rcParams['xtick.minor.size'] = 2
rcParams['ytick.major.size'] = 4
rcParams['ytick.minor.size'] = 2

au = 1.5e13
msun = 2.0e33
G = 6.67e-8
m_star = 31.7 * msun

# Mass of the perturber
m_sl = 3.2 * msun

def plot_one_disk(r_encounter, v_encounter, ax1, title=''):
"""
Encounter parameters based on the formalism by
D’Onghia, E., Vogelsberger, M., Faucher-Giguere, C.-A. & Hernquist, L. 
"Quasi-resonant Theory of Tidal Interactions".
2010 ApJ,725,353
"""
    b = r_encounter * au   #
    rout = 2000 * au
    v_sl = v_encounter * 1e5
    rr = np.logspace(np.log10(rout) - 3, np.log10(rout), 1000)

    omega = (G * m_star / rr**3)**0.5
    alpha = omega * b / v_sl

    term_o1_p = 2 * G * m_sl / b**2 / v_sl * rr
    term_o1_m = (2 * np.pi)**0.5 * alpha**(1.5) * np.exp(- alpha)
    term_o1 = term_o1_p * term_o1_m

    term_o2_p = G * m_sl / b**3 / v_sl * rr ** 2
    term_o2_m = 2 * alpha**2 * (1 + 4 * alpha) * kn(0, 2 * alpha)
    term_o2 = term_o2_p * term_o2_m

    ax1.plot(rr / au, term_o1 / 1e5, 'r--', label='1st-order term')
    ax1.plot(rr / au, term_o2 / 1e5, 'b--', label='2nd-order term')
    ax1.plot(rr / au, (term_o1 + term_o2) / 1e5, 'k-', label='1st+2nd-order')
    ax1.set_xlim(0, 2000.)
    ax1.set_ylim(0, 1)
    ax1.set_title(title)

if __name__ == '__main__':

    # Angular velocity of the perturber at the periastron (rad s-1)
    oo = np.array([1.0e-11,1.3e-11,1.6e-11])
    xrange = len(oo)

    figure1 = plt.figure(figsize=(6,6))
    outer = gs.GridSpec(1,1)
    outer.update(left=0.16,right=0.98,top=0.90,bottom=0.11)

    # Periastron distance of the perturber (AU)
    pp = [2500, 2000, 1500]

    count = 1
    for countj, j in enumerate(pp):
        for counti, i in enumerate(range(xrange)):
            ax1 = figure1.add_subplot(len(pp), xrange, count)
            ax1.set_xlabel('$r$ (AU)')
            ax1.set_xticks([0,500,1000,1500,2000])
            ax1.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
            if countj < len(pp) - 1:
                ax1.set_xlabel(' ')
                ax1.set_xticklabels(' ')
            if counti > 0:
                ax1.set_ylabel(' ')
                ax1.set_yticklabels(' ')
            r = j
            omega = oo[i]
            v = omega * r * au / 1e5
            if countj == 0:
                if counti == 0:
                    omega_s = r'$\omega_\mathrm{peri}$=' + "{:.1f}".format(omega*1e11) + r'$\times$$10^{-11}$$\;\rm s^{-1}$'
                else:
                    omega_s = "{:.1f}".format(omega*1e11) + r'$\times$$10^{-11}$$\;\rm s^{-1}$'
            else:
                omega_s = ' '

            v_0  = (G * m_star / r / au) ** 0.5
            omega_0 =  v_0 / r / au
            
            plot_one_disk(r, v, ax1, omega_s)            
            ax1.text(0.05, 0.1, r"$\omega_\mathrm{peri} / \omega_{\rm Kep}$ = " +
                     "{:.2f}".format(omega/omega_0), transform=ax1.transAxes)
            
            ax1a = ax1.twinx()
            ax1a.set_yticks([])
            ax1.tick_params(left=True, right=True, bottom=True, top=True)

            if counti == xrange -1:
                if countj == len(pp)-1:
                    ax1a.set_ylabel(r'$r_\mathrm{peri}$=' + str(j) + ' AU', fontsize=12)
                else:
                    ax1a.set_ylabel(str(j) + ' AU', fontsize=12)

            if counti == 0:
                ax1.set_ylabel(r'$\delta_{\rm v}\rm (km/s)$')
                if count == 1:
                    ax1.legend(loc='upper left', fontsize=9)
            count = count + 1
    figure1.subplots_adjust(left=0.08,right=0.95,bottom=0.07,top=0.95)
    figure1.savefig('perturbed_disk_3x3.pdf',dpi=300)

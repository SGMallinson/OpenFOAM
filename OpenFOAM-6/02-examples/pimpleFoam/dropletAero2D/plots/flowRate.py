#!/usr/bin/env python

"""Plot outlet flowRate
versus time and compare with Couette value for MSR-368

:author: S.G. Mallinson <sam>

:created: 2019-10-18

"""

from matplotlib.pyplot import subplots
import matplotlib.ticker as ticker
from pandas import read_csv
from numpy import loadtxt, min, max


def main(name: str,
         outlet: str) -> None:

    Up = 350 # mm/s
    PPS = 1.5 # mm
    wdom = 10.0 * PPS
    Qc = 0.5 * wdom * PPS * Up # uL / s
    nu = 15.5 # mm^2 / s
    tEst = PPS**2 / nu
    print(Qc)
    tOutlet, Qoutlet = loadtxt(outlet, skiprows=4, usecols=(0, 1),
                               dtype = 'float32', unpack='True')

    fig, ax = subplots(2)

    ax[0].plot(tOutlet, 1e9 * Qoutlet,
            color='red', linestyle='-', marker='None',
            label='outlet')
    #ax[0].plot([0, max(tOutlet)], [Qc, Qc],
    ax[0].plot([0, 2], [Qc, Qc],
            color='black', linestyle='-', marker='None',
            label='Couette')
    ax[0].plot([tEst, tEst], [0, Qc],
            color='green', linestyle='-', marker='None',
            label=r'$h^2 / \nu$')
    ax[0].set_ylabel(r'$Q$, uL / s')
    ax[0].set_xlim([0, 0.2])
    #ax[0].set_ylim([0, 4000])
    ax[0].legend(loc='best')

    ax[1].semilogy(tOutlet, 1.0 - Qoutlet / Qoutlet[-1],
            color='blue', linestyle='-', marker='None',
            label='inlet')
    ax[1].set_ylabel(r'shortfall from final value')
    ax[1].set_xlabel('$t$, s')
    ax[1].set_xlim([0, 0.2])

    fig.savefig(name + '.png')
    

if __name__ == '__main__':

    from sys import argv
    from os.path import basename, splitext

    main(splitext(basename(argv[0]))[0], *argv[1:])

#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

import argparse
import os

import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op


def main():
    options = _parse_args()

    for run in options.run:
        shortname = os.path.basename(os.path.abspath(run))

        energy = np.atleast_2d(np.loadtxt(os.path.join(run, 'energy.tsv')))
        plaquette = np.atleast_2d(np.loadtxt(os.path.join(run, 'plaquette.tsv')))
        energy_reject = np.atleast_2d(np.loadtxt(os.path.join(run, 'energy-reject.tsv')))
        plaquette_reject = np.atleast_2d(np.loadtxt(os.path.join(run, 'plaquette-reject.tsv')))

        boltzmann = np.loadtxt(os.path.join(run, 'boltzmann.tsv'))

        fig = pl.figure()

        ax_energy = fig.add_subplot(2, 1, 1)
        if len(energy) > 1:
            ax_energy.plot(energy[:, 0], energy[:, 1], color='green')
        if len(energy_reject) > 2:
            ax_energy.plot(energy_reject[:, 0], energy_reject[:, 1],
                           linestyle='none', color='red', marker='.', alpha=0.3)
        ax_energy.set_title(shortname)
        ax_energy.set_xlabel('Markov Chain Index')
        ax_energy.set_ylabel(r'Average Energy')
        ax_energy.grid(True)
        ax_energy.margins(0.05)

        ax_plaquette = fig.add_subplot(2, 1, 2)
        if len(plaquette) > 1:
            ax_plaquette.plot(plaquette[:, 0], plaquette[:, 1], color='green')
        if len(plaquette_reject) > 2:
            ax_plaquette.plot(plaquette_reject[:-1, 0], plaquette_reject[:-1, 1],
                              linestyle='none', color='red', marker='.', alpha=0.3)
        ax_plaquette.set_xlabel('Markov Chain Index')
        ax_plaquette.set_ylabel(r'Average Plaquette')
        ax_plaquette.grid(True)
        ax_plaquette.margins(0.05)

        fig.tight_layout()
        fig.savefig('plot-averages-{}.pdf'.format(shortname))

        fig = pl.figure()
        ax = fig.add_subplot(2, 1, 1)
        ax.hist(boltzmann, bins=30)
        ax.set_title(shortname)
        ax.set_xlabel('$\Delta H$')
        ax.set_ylabel('Frequency')
        ax.grid(True)
        ax = fig.add_subplot(2, 1, 2)
        ax.plot(boltzmann)
        ax.set_xlabel('Tried trajectory')
        ax.set_ylabel('$\Delta H$')
        ax.grid(True)
        fig.tight_layout()
        fig.savefig('plot-energy_diff-{}.pdf'.format(shortname))

        print('Mean of exp( ΔH):', np.mean(np.exp(boltzmann)))
        print('Mean of exp(−ΔH):', np.mean(np.exp(-boltzmann)))


def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :r
    type: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('run', nargs='+')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()

#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

import argparse
import configparser
import os

import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op


def linear(x, a, b):
    return a * x + b


def main():
    options = _parse_args()

    fig = pl.figure()
    ax_plaquette = fig.add_subplot(1, 2, 1)
    ax_beta = fig.add_subplot(1, 2, 2)

    beta_x = []
    beta_y = []
    beta_y_err = []

    for run in options.run:
        config = configparser.ConfigParser()
        config.read(os.path.join(run, 'hmc.ini'))
        beta = float(config['md']['beta'])

        data = np.atleast_2d(np.loadtxt(os.path.join(run, 'plaquette.tsv')))
        if len(data) > 1:
            ax_plaquette.plot(data[:, 0], data[:, 1], label=r'$\beta$ = {}'.format(beta))

            if len(data[:, 1]) > 20:
                beta_x.append(beta)
                beta_y.append(np.mean(data[20:, 1]))
                beta_y_err.append(np.std(data[20:, 1]))

    beta_x = np.array(beta_x)
    beta_y = np.array(beta_y)
    beta_y_err = np.array(beta_y_err)

    ax_plaquette.set_title('$8^4$ lattice, warm start\nstep 0.01, steps 100')
    ax_plaquette.grid(True)
    ax_plaquette.legend(loc='best', prop={'size': 8})
    ax_plaquette.margins(0.05)
    ax_plaquette.set_xlabel(r'Markov Time $\tau_\mathrm{M}$')
    ax_plaquette.set_ylabel('W(1×1)')

    ax_beta.errorbar(beta_x, beta_y, yerr=beta_y_err, marker='o', linestyle='none')

    sel = beta_x <= 4
    popt, pconv = op.curve_fit(linear, beta_x[sel], beta_y[sel], sigma=beta_y_err[sel])
    x = np.linspace(0, 2, 2)
    y = linear(x, *popt)
    ax_beta.plot(x, y)

    np.savetxt('plaquette-vs-beta.tsv',
               np.column_stack([beta_x, beta_y, beta_y_err]))

    print(*popt)

    ax_beta.grid(True)
    ax_beta.margins(0.05)
    ax_beta.set_xlabel(r'$\beta$')
    ax_beta.set_ylabel(r'W(1×1) average over $20 \leq \tau_\mathrm{M}$')

    fig.tight_layout()
    fig.savefig('plot-beta.pdf')


def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :rtype: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('run', nargs='+')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()

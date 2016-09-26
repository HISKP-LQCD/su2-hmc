#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

import argparse
import configparser
import os

import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op

import bootstrap

def main():
    options = _parse_args()

    steps_list = []
    rate_val_list = []
    rate_err_list = []

    for run in options.run:
        config = configparser.ConfigParser()
        config.read(os.path.join(run, 'hmc.ini'))

        beta = float(config['md']['beta'])
        steps = float(config['md']['time_step'])

        acceptance = np.loadtxt(os.path.join(run, 'accept.tsv'))

        start = 30
        values = acceptance[start:]

        val, err = bootstrap.bootstrap_and_transform(np.mean, values)

        print(steps, val, err)

        steps_list.append(steps)
        rate_val_list.append(val)
        rate_err_list.append(err)

    pl.errorbar(steps_list, rate_val_list, yerr=rate_err_list, marker='+', linestyle='none')
    pl.xscale('log')
    pl.xlabel(r'Time Step $\Delta \tau$ in Molecular Dynamics')
    pl.ylabel('Acceptance Rate')
    pl.grid(True)
    pl.margins(0.05)
    pl.tight_layout()
    pl.savefig('acceptance-rate.pdf')

    np.savetxt('acceptance-vs-time-step.tsv',
               np.column_stack([steps_list, rate_val_list, rate_err_list]))




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

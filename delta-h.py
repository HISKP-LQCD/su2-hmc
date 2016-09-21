#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

import argparse
import random

import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op
import scipy as sp


def autocorrelation_time(data):
    length = len(data)
    correlated = np.array([
        np.mean(data * np.roll(data, i))
        for i in range(length // 2)
    ]) - np.mean(data)**2
    autocorrelation = correlated / correlated[0]

    fig = pl.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('$t$')
    ax.set_ylabel(r'Autocorrelation $\Gamma(t)$')
    ax.plot(autocorrelation)
    ax.margins(0.05)
    fig.tight_layout()
    fig.savefig('corr.pdf')

    integral = sp.integrate.trapz(autocorrelation)
    print(integral)
    integral = sp.integrate.quad(autocorrelation)
    print(integral)
    print(np.sum(autocorrelation))


def main():
    options = _parse_args()

    delta_h = np.loadtxt(options.filename)

    boltzmann = np.exp(- delta_h)[100:]

    autocorrelation_time(boltzmann)

    print(np.mean(boltzmann), np.std(boltzmann))

    block_lens = list(range(1, 101))
    mm = []
    ms = []
    sm = []
    ss = []

    for block_len in block_lens:
        divisible = len(boltzmann) // block_len * block_len
        reshaped = boltzmann[:divisible].reshape((-1, block_len))
        blocked = np.mean(reshaped, axis=1)
        means = []
        stds = []
        for sample in range(3 * len(blocked)):
            sampled = [random.choice(blocked) for i in range(len(blocked))]
            means.append(np.mean(sampled))
            stds.append(np.std(sampled))
        mm.append(np.mean(means))
        sm.append(np.std(means))
        ms.append(np.mean(stds))
        ss.append(np.std(stds))

    fig = pl.figure()

    ax = fig.add_subplot(2, 1, 1)
    ax.errorbar(block_lens, mm, yerr=sm)
    ax.set_xlabel('Block size')
    ax.set_ylabel('Mean of blocked mean$')
    ax.margins(0.05)

    ax = fig.add_subplot(2, 1, 2)
    ax.errorbar(block_lens, ms, yerr=ss)
    ax.set_xlabel('Block size')
    ax.set_ylabel('Mean of blocked std')
    ax.margins(0.05)

    fig.tight_layout()
    fig.savefig('1-errs.pdf')

    fig = pl.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(boltzmann)
    fig.savefig('delta-h.pdf')


def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :rtype: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filename')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()

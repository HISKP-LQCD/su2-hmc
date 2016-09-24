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
    ax.grid(True)
    fig.tight_layout()
    fig.savefig('corr.pdf')

    integral = sp.integrate.trapz(autocorrelation)
    print(integral)
    print(np.sum(autocorrelation))


def main():
    options = _parse_args()

    data = np.loadtxt(options.filename)
    shape = data.shape
    if len(shape) == 1:
        delta_h = data
    else:
        delta_h = data[:, 1]

    boltzmann = np.exp(- delta_h)[20:]

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

    mm = np.array(mm)
    ms = np.array(ms)
    sm = np.array(sm)
    ss = np.array(ss)

    fig = pl.figure()

    ax = fig.add_subplot(2, 2, 1)
    #ax.errorbar(block_lens, mm, yerr=sm)
    ax.fill_between(block_lens, mm - sm, mm + sm, alpha=0.3, linewidth=0)
    ax.plot(block_lens, mm)
    ax.set_xlabel('Block size')
    ax.set_ylabel('Mean of blocked mean')
    ax.margins(0.05)
    ax.grid(True)

    ax = fig.add_subplot(2, 2, 3)
    ax.plot(block_lens, sm)
    ax.set_xlabel('Block size')
    ax.set_ylabel('Std of blocked mean')
    ax.margins(0.05)
    ax.grid(True)

    ax = fig.add_subplot(2, 2, 2)
    #ax.errorbar(block_lens, ms, yerr=ss)
    ax.fill_between(block_lens, ms - ss, ms + ss, alpha=0.3, linewidth=0)
    ax.plot(block_lens, ms)
    ax.set_xlabel('Block size')
    ax.set_ylabel('Mean of blocked std')
    ax.margins(0.05)
    ax.grid(True)

    ax = fig.add_subplot(2, 2, 4)
    ax.plot(block_lens, ss)
    ax.set_xlabel('Block size')
    ax.set_ylabel('Std of blocked std')
    ax.margins(0.05)
    ax.grid(True)

    fig.tight_layout()
    fig.savefig('1-errs.pdf')

    fig = pl.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(boltzmann)
    ax.set_xlabel(r'$\exp(-\Delta H)$')
    ax.set_ylabel('Counts')
    fig.tight_layout()
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

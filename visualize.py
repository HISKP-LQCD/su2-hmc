#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

import argparse

import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D



def main():
    options = _parse_args()

    raw = np.loadtxt(options.numbers)
    data = raw.reshape((8, 8, 8, 8))

    print(np.min(data), np.max(data))

    verts, faces = measure.marching_cubes(data[0, :, :, :], 0)

    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2], cmap='Spectral', lw=1)
    fig.savefig('test.png')
    
def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :rtype: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('numbers')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()

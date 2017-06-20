#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

import argparse
import os

import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
import jinja2

template_text = r'''

#include "colors.inc"

#declare damped = 0.4;

camera {
    location <-5, -5, -1>
    look_at <4, 4, 4>
}

global_settings {
    ambient_light Gray10
}

light_source {
    <-2, 0, 0>
    color rgb<damped, damped, 1.0>
    area_light <0, 8, 0> <0, 0, 8> 3 3
}

light_source {
    <0, -2, 0>
    color rgb<damped, 1, damped>
    area_light <8, 0, 0> <0, 0, 8> 3 3
}

light_source {
    <0, 0, -2>
    color rgb<1, damped, damped>
    area_light <0, 8, 0> <8, 0, 0> 3 3
}

blob {
    threshold {{ threshold }}

    {% for i in i_list -%}
    {% for j in j_list -%}
    {% for k in k_list -%}
    sphere {
        <{{ (i+L)%L }}, {{ (j+L)%L }}, {{ (k+L)%L }}>,
        1,
        {{ data[(i+L)%L, (j+L)%L, (k+L)%L] }}
    }
    {% endfor -%}
    {% endfor -%}
    {% endfor -%}

    texture{
        pigment{
            color rgb<1.0, 1.0, 1.0>
        }
        finish {
            phong 1
        }
    }
    
    no_shadow
}

'''

def main():
    options = _parse_args()

    raw = np.loadtxt(options.numbers)
    L = int(len(raw)**(1/4))
    data = raw.reshape((L, L, L, L))  

    print(np.min(data), np.max(data))

    pl.hist(raw)
    pl.savefig('hist.pdf')

    template = jinja2.Template(template_text)

    shape = data.shape

    base, ext = os.path.splitext(options.numbers)
    basename = os.path.basename(base)

    for t in range(shape[0]):
        use = data[t, :, :, :]
        rendered = template.render(
            i_list=list(range(shape[1])),
            j_list=list(range(shape[2])),
            k_list=list(range(shape[3])),
            data=use,
            threshold=np.percentile(data, 99),
            L=L,
        )

        with open('scene-{}-t{:03d}.pov'.format(basename, t), 'w') as f:
            f.write(rendered)
    
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

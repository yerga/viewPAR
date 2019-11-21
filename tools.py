# -*- coding: utf-8 -*-
import numpy as np


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]


def convert_reference(potentials, usedref, toref, ph):
    # RefElectrodes:  0:RHE, 1:AgAgCl

    if toref == 0 and usedref == 1:
        for i in range(len(potentials)):
            potentials[i] = np.array(potentials[i]).astype(float) + 0.059*ph + 0.197

    return potentials

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# This module contains code to plot
# PBF scode file
##################################

# standard library
import sys

# 3rd party library
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

sns.set()

ZERO_POWER = 1.0

def plot_roads_3D(roads):
    """plot roads in 3D."""
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    for x1, y1, x2, y2, z, power, speed in roads:
        # ax.plot([x1, x2], [y1, y2], [z, z], 'b-')
        if power > ZERO_POWER:
            ax.plot([x1, x2], [y1, y2], [z, z], 'r-')
        else:
            ax.plot([x1, x2], [y1, y2], [z, z], 'b-')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    return fig, ax


def plot_roads_2D(roads):
    """plot roads in 2D."""
    fig, ax = plt.subplots(figsize=(8, 8))
    for x1, y1, x2, y2, z, power, speed in roads:
        if power > ZERO_POWER:
            ax.plot([x1, x2], [y1, y2], 'r-')
        else:
            ax.plot([x1, x2], [y1, y2], 'b-')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    return fig, ax


if __name__ == "__main__":
    if len(sys.argv) != 2 or (not sys.argv[1].endswith('.scode')):
        print("Usage: >> python {} {}".format(sys.argv[0], "<scode file>"))
        sys.exit(1)
    path_file = sys.argv[1]
    roads = np.loadtxt(path_file)
    nlayers = len(set([z for z in roads[:, 4]]))
    if nlayers == 1:
        fig, ax = plot_roads_2D(roads)
    else:
        fig, ax = plot_roads_3D(roads)
    plt.show()

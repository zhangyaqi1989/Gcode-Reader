#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# refactor the original version of
# Gcode reader
##################################

# standard library
import argparse
import collections
from enum import Enum
import math
import os.path
import pprint
import statistics
import sys

# 3rd party library
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import seaborn as sns


class LayerValueError(Exception):
    """layer Value error."""
    pass


class GcodeType(Enum):
    """ enum of GcodeType """

    FDM_REGULAR = 1
    FDM_STRATASYS = 2
    LPBF_REGULAR = 3
    LPBF_SCODE = 4


def save_figure(fig, filename, dpi=100):
    """save figure to a file."""
    _, ext = filename.rsplit('.', 1)
    fig.savefig(filename, format=ext, dpi=dpi, bbox_inches='tight')
    print('Saving to {:s} with {:d} DPI'.format(filename, dpi))


def create_fig_axis(figsize=(8, 8), projection='2d'):
    """create and return fig, ax based on figure size and projection."""
    if projection not in ['2d', '3d', '2D', '3D']:
        raise ValueError(f"invalid projection: {projection}.")
    projection = projection.lower()
    if projection == '2d':
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
    return fig, ax


def create_movie_writer(title='Movie Writer', fps=15):
    """return a ffmpeg writer."""
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title=title, artist='Matplotlib',
            comment='Movie Support')
    writer = FFMpegWriter(fps=fps, metadata=metadata)
    return writer


def get_parser():
    """set up parser and return it."""
    parser = argparse.ArgumentParser(description='G-Code Reader')
    parser.add_argument(dest='gcode_file', action='store',
            help='path of the input G-Code file')
    parser.add_argument('-t', '--type', dest='file_type', type=int, default=1, choices=[1,2,3,4],
            help='1: Regular FDM; 2: Stratasys FDM;\n 3: Regular LPBF; 4: Scode LPBF')
    parser.add_argument('-lo', '--low', dest='low_layer', type=int, default=1,
            help='lowest layer number')
    parser.add_argument('-hi', '--high', dest='high_layer', type=int, default=10000,
            help='highest layer number')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-p", "--plot", action="store_true", help='plot the part')
    group.add_argument("-a", "--animate", action="store_true", help='animate the printing process')
    group.add_argument("-m", "--mesh", action="store_true", help='plot the meshing of the path')
    parser.add_argument('-s', '--save', dest='out_file', type=str,
            action='store', help='path of output file')
    return parser


def command_line_runner():
    """command line runner."""
    # 1. parse arguments
    parser = get_parser()
    args = parser.parse_args()

    print(type(args.file_type))
    print(type(args.low_layer), type(args.high_layer))
    print(args.gcode_file)
    print(args.plot, args.animate, args.mesh)
    print(args.out_file) # if not set, args.out_file is None

    # 2. run code based on args
    gcode_reader = GcodeReader(args.gcode_file, args.file_type)
    if args.plot:
        fig, ax = gcode_reader.plot(low_layer=args.low_layer, high_layer=args.high_layer)
    if args.mesh:
        fig, ax = gcode_reader.mesh(low_layer=args.low_layer, high_layer=args.high_layer)
    if args.animate:
        fig, ax = gcode_reader.animate(low_layer=args.low_layer, high_layer=args.high_layer)


class GcodeReader:
    """GCode Reader."""

    def __init__(self, gcode_file, file_type):
        self.gcode_file = gcode_file
        self.file_type = GcodeType(file_type)
        self._read_gcode()

    def _read_gcode(self):
        """read G-Code file based on file_type."""
        if self.file_type == GcodeType.FDM_REGULAR:
            self._read_fdm_regular()
        elif self.file_type == GcodeType.FDM_STRATASYS:
            self._read_fdm_stratasys()
        elif self.file_type == GcodeType.LPBF_REGULAR:
            self._read_lpbf_regular()
        elif self.file_type == GcodeType.LPBF_SCODE:
            self._read_lpbf_scode()
        else:
            assert False, f"Unsupport file type {self.file_type}."

    def _read_fdm_regular(self):
        pass

    def _read_fdm_stratasys(self):
        pass

    def _read_lpbf_regular(self):
        pass

    def _read_lpbf_scode(self):
        pass

    def plot(self, low_layer, high_layer):
        """plot part from low_layer to high_layer."""
        print(f"Plot part from {low_layer} to {high_layer}.")
        return None, None

    def mesh(self, low_layer, high_layer):
        """mesh the part and plot the mesh from low_layer to high_layer."""
        print(f"Plot meshing of part from {low_layer} to {high_layer}.")
        return None, None

    def animate(self, low_layer, high_layer):
        """animate the printing of part from low_layer to high_layer."""
        print(f"Animate printing of part from {low_layer} to {high_layer}.")
        return None, None


if __name__ == "__main__":
    command_line_runner()
    # fig, ax = create_fig_axis((8, 8), '2d')
    # writer = create_movie_writer()

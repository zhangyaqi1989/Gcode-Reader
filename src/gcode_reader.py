#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# Gcode reader
##################################
# TODO:
# 1. number of subpaths
# 2. number of nozzle travels
# 3. total distance
# 4. number of elements in each layer
# 5. support stratasys file
# 6. support LPBF file
##################################

# standard library
import math
import os.path
import sys
import argparse
from enum import Enum
import pprint

# third party library
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

pp = pprint.PrettyPrinter(indent=4)


class GcodeType(Enum):
    """ enum of GcodeType """

    FDM_REGULAR = 1
    FDM_STRATASYS = 2
    LPBF = 3

    @classmethod
    def has_value(cls, value):
        return any(value == item.value for item in cls)


class GcodeReader:
    """ Gcode reader class """

    def __init__(self, filename, filetype=GcodeType.FDM_REGULAR):
        if not os.path.exists(filename):
            print("{} does not exist!".format(filename))
            sys.exit(1)
        self.filename = filename
        self.filetype = filetype
        # print(self.filetype)
        self.n_segs = 0
        self.segs = None
        self.n_layers = 0
        self.index_bars = []
        self.summary = None
        self.lengths = None
        self.subpaths = None
        # read file to populate variables
        self._read()

    def _read(self):
        """ read the file and stores segs into data structure """
        if self.filetype == GcodeType.FDM_REGULAR:
            self._read_fdm_regular()
        else:
            print("file type is not supported")
            sys.exit(1)

    def _read_fdm_regular(self):
        """ read fDM regular gcode type """
        with open(self.filename) as infile:
            # read nonempty lines
            lines = [line.strip() for line in infile.readlines()
                     if line.strip()]
            # only keep line that starts with 'G1'
            lines = [line for line in lines if line.startswith('G1')]
        # pp.pprint(lines) # for debug
        self.segs = []
        temp = -float('inf')
        gxyzef = [temp, temp, temp, temp, temp, temp]
        d = dict(zip(['G', 'X', 'Y', 'Z', 'E', 'F'], range(6)))
        for i, line in enumerate(lines):
            old_gxyzef = gxyzef[:]
            for token in line.split():
                gxyzef[d[token[0]]] = float(token[1:])
            if gxyzef[3] > old_gxyzef[3]:
                self.n_layers += 1
                self.index_bars.append(i)
            if (gxyzef[0] == 1 and gxyzef[1:3] != old_gxyzef[1:3]
                    and gxyzef[3] == old_gxyzef[3]
                    and gxyzef[4] > old_gxyzef[4]):
                x0, y0, z = old_gxyzef[1:4]
                x1, y1, _ = gxyzef[1:4]
                self.segs.append((x0, y0, x1, y1, z))
        self.n_segs = len(self.segs)
        self.segs = np.array(self.segs)
        self.index_bars.append(self.n_segs)
        assert(len(self.index_bars) - self.n_layers == 1)

    def _compute_subpaths(self):
        """ compute subpaths
            a subpath is represented by (xs, ys, zs)
        """
        if not self.subpaths:
            self.subpaths = []
            x0, y0, x1, y1, z = self.segs[0, :]
            xs, ys, zs = [x0, x1], [y0, y1], [z, z]
            for x0, y0, x1, y1, z in self.segs[1:, :]:
                if x0 != xs[-1] or y0 != ys[-1] or z != zs[-1]:
                    self.subpaths.append((xs, ys, zs))
                    xs, ys, zs = [x0, x1], [y0, y1], [z, z]
                else:
                    xs.append(x1)
                    ys.append(y1)
                    zs.append(z)
            if len(xs) != 0:
                self.subpaths.append((xs, ys, zs))

    def plot(self, color='blue'):
        """ plot the whole part in 3D """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        assert(self.n_segs > 0)
        self._compute_subpaths()
        for xs, ys, zs in self.subpaths:
            ax.plot(xs, ys, zs)
        plt.show()

    def plot_layer(self, layer=1):
        """ plot a specific layer in 2D """
        # make sure layer is in [1, self.n_layers]
        layer = max(layer, 1)
        layer = min(self.n_layers, layer)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        self._compute_subpaths()
        left, right = self.index_bars[layer - 1], self.index_bars[layer]
        for xs, ys, _ in self.subpaths[left: right]:
            ax.plot(xs, ys)
        plt.show()

    def describe(self):
        if not self.summary:
            self.lengths = [math.hypot(x1 - x0, y1 - y0) for x0, y0, x1, y1, _
                            in self.segs]
            series = pd.Series(self.lengths)
            self.summary = series.describe()
        print("1. Line segments information: ")
        print(self.summary)


def command_line_runner():
    """ main function """
    # parse arguments
    parser = argparse.ArgumentParser(description='Gcode Reader')
    parser.add_argument(dest='gcode_file', help='gcode file', action='store')
    parser.add_argument('-t', dest='filetype', help="""File Type
            1: regular FDM; 2: Stratasys FDM; 3: LPBF""",
                        required=True, type=int, action='store')
    args = parser.parse_args()
    # handle Gcode file type
    if not GcodeType.has_value(args.filetype):
        print('Invalid G-code file type: {:d}'.format(args.filetype))
        print('Valid types are listed below')
        for gcode_type in GcodeType:
            print('{:s} : {:d}'.format(gcode_type.name, gcode_type.value))
        sys.exit(1)
    else:
        filetype = GcodeType(args.filetype)
    gcode_reader = GcodeReader(filename=args.gcode_file, filetype=filetype)
    gcode_reader.describe()
    gcode_reader.plot_layer(layer=1)
    # gcode_reader.plot()


if __name__ == "__main__":
    command_line_runner()

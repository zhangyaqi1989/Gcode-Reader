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


class LayerError(Exception):
    """ layer number error """
    pass


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
        self.seg_index_bars = []
        self.subpath_index_bars = []
        self.summary = None
        self.lengths = None
        self.subpaths = None
        # read file to populate variables
        self._read()

    def _read(self):
        """ read the file and stores segs into data structure """
        if self.filetype == GcodeType.FDM_REGULAR:
            self._read_fdm_regular()
        elif self.filetype == GcodeType.FDM_STRATASYS:
            self._read_fdm_stratasys()
        elif self.filetype == GcodeType.LPBF:
            self._read_lpbf()
        else:
            print("file type is not supported")
            sys.exit(1)

    def _read_lpbf(self):
        """ read LPBF gcode """
        with open(self.filename) as infile:
            # read nonempty lines
            lines = [line.strip() for line in infile.readlines()
                     if line.strip()]
            # only keep line that starts with 'N'
            lines = [line for line in lines if line.startswith('N')]
        # pp.pprint(lines) # for debug
        self.segs = []
        self.powers = []
        temp = -float('inf')
        ngxyzfl = [temp, temp, temp, temp, temp, temp, temp]
        d = dict(zip(['N', 'G', 'X', 'Y', 'Z', 'F', 'L'], range(7)))
        seg_count = 0
        for line in lines:
            old_ngxyzfl = ngxyzfl[:]
            tokens = line.split()
            for token in tokens:
                ngxyzfl[d[token[0]]] = float(token[1:])
            if ngxyzfl[d['Z']] > old_ngxyzfl[d['Z']]:
                self.n_layers += 1
                self.seg_index_bars.append(seg_count)
            if (ngxyzfl[1] == 1 and ngxyzfl[2:4] != old_ngxyzfl[2:4]
                    and ngxyzfl[4] == old_ngxyzfl[4]
                    and ngxyzfl[5] > 0):
                x0, y0, z = old_ngxyzfl[2:5]
                x1, y1 = ngxyzfl[2:4]
                self.segs.append((x0, y0, x1, y1, z))
                seg_count += 1
        self.n_segs = len(self.segs)
        self.segs = np.array(self.segs)
        self.seg_index_bars.append(self.n_segs)
        # print(self.n_layers)
        assert(len(self.seg_index_bars) - self.n_layers == 1)

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
        seg_count = 0
        for line in lines:
            old_gxyzef = gxyzef[:]
            for token in line.split():
                gxyzef[d[token[0]]] = float(token[1:])
            if gxyzef[3] > old_gxyzef[3]:  # z value
                self.n_layers += 1
                self.seg_index_bars.append(seg_count)
            if (gxyzef[0] == 1 and gxyzef[1:3] != old_gxyzef[1:3]
                    and gxyzef[3] == old_gxyzef[3]
                    and gxyzef[4] > old_gxyzef[4]):
                x0, y0, z = old_gxyzef[1:4]
                x1, y1 = gxyzef[1:3]
                self.segs.append((x0, y0, x1, y1, z))
                seg_count += 1
        self.n_segs = len(self.segs)
        self.segs = np.array(self.segs)
        self.seg_index_bars.append(self.n_segs)
        assert(len(self.seg_index_bars) - self.n_layers == 1)

    def _read_fdm_stratasys(self):
        """ read stratasys fdm G-code file """
        self.areas = []
        self.is_supports = []
        self.styles = []
        self.deltTs = []
        self.segs = []
        temp = -float('inf')
        # x, y, z, area, deltaT, is_support, style
        xyzATPS = [temp, temp, temp, temp, temp, False, '']
        seg_count = 0
        with open(self.filename, 'r') as in_file:
            lines = in_file.readlines()
            # means position denoted by the line is the start of subpath
            is_start = True
            for line in lines:
                if line.startswith('#'):
                    continue
                if not line.strip():  # skip empty line
                    start = True
                    continue
                old_xyzATPS = xyzATPS[:]
                tokens = line.split()
                # print(tokens)
                xyzATPS[:5] = [float(token) for token in tokens[:5]]
                xyzATPS[5] = bool(tokens[5])
                xyzATPS[6] = tokens[6]
                if xyzATPS[2] != old_xyzATPS[2]:  # z value
                    self.seg_index_bars.append(seg_count)
                    self.n_layers += 1
                elif not start:
                    # make sure is_support and style do not change
                    assert(xyzATPS[5:] == old_xyzATPS[5:])
                    x0, y0 = old_xyzATPS[:2]
                    x1, y1, z = xyzATPS[:3]
                    self.segs.append((x0, y0, x1, y1, z))
                    seg_count += 1
                    self.areas.append(xyzATPS[3])
                    self.deltTs.append(xyzATPS[4])
                    self.is_supports.append(xyzATPS[5])
                    self.styles.append(xyzATPS[6])
                start = False
            self.n_segs = len(self.segs)
            self.segs = np.array(self.segs)
            self.seg_index_bars.append(self.n_segs)
            # print(self.seg_index_bars)

    def _compute_subpaths(self):
        """ compute subpaths
            a subpath is represented by (xs, ys, zs)
        """
        if not self.subpaths:
            self.subpaths = []
            self.subpath_index_bars = [0]
            x0, y0, x1, y1, z = self.segs[0, :]
            xs, ys, zs = [x0, x1], [y0, y1], [z, z]
            for x0, y0, x1, y1, z in self.segs[1:, :]:
                if x0 != xs[-1] or y0 != ys[-1] or z != zs[-1]:
                    self.subpaths.append((xs, ys, zs))
                    if z != zs[-1]:
                        self.subpath_index_bars.append(len(self.subpaths))
                    xs, ys, zs = [x0, x1], [y0, y1], [z, z]
                else:
                    xs.append(x1)
                    ys.append(y1)
                    zs.append(z)
            if len(xs) != 0:
                self.subpaths.append((xs, ys, zs))
            self.subpath_index_bars.append(len(self.subpaths))
            # print(self.subpath_index_bars)
            # print(self.segs)

    def plot(self, color='blue'):
        """ plot the whole part in 3D """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        assert(self.n_segs > 0)
        self._compute_subpaths()
        for xs, ys, zs in self.subpaths:
            ax.plot(xs, ys, zs)
        plt.show()

    def plot_layers(self, min_layer, max_layer):
        """ plot the layers in [min_layer, max_layer) in 3D """
        if (min_layer >= max_layer or min_layer < 1 or max_layer >
                self.n_layers):
            raise LayerError("Layer number is invalid!")
        self._compute_subpaths()
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        left, right = (self.subpath_index_bars[min_layer - 1],
                       self.subpath_index_bars[max_layer - 1])
        for xs, ys, zs in self.subpaths[left: right]:
            ax.plot(xs, ys, zs)
        plt.show()

    def plot_layer(self, layer=1):
        """ plot a specific layer in 2D """
        # make sure layer is in [1, self.n_layers]
        # layer = max(layer, 1)
        # layer = min(self.n_layers, layer)
        if layer < 1 or layer > self.n_layers:
            raise LayerError("Layer number is invalid!")
        self._compute_subpaths()
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        left, right = (self.subpath_index_bars[layer - 1],
                       self.subpath_index_bars[layer])
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
        print("2. number of layers: {:d}".format(self.n_layers))
        self._compute_subpaths()
        # print(len(self.seg_index_bars))
        # print(len(self.subpath_index_bars))
        data = {'# segments': np.array(self.seg_index_bars[1:]) -
                np.array(self.seg_index_bars[:-1]),
                'layer': np.arange(1, self.n_layers + 1),
                '# subpaths': np.array(self.subpath_index_bars[1:]) -
                np.array(self.subpath_index_bars[:-1]),
                }
        df = pd.DataFrame(data)
        df = df.set_index('layer')
        print(df)


def command_line_runner():
    """ main function """
    # parse arguments
    parser = argparse.ArgumentParser(description='Gcode Reader')
    parser.add_argument(dest='gcode_file', help='gcode file', action='store')
    parser.add_argument('-t', '--type', dest='filetype', help="""File Type
            1: regular FDM; 2: Stratasys FDM; 3: LPBF""",
                        required=True, type=int, action='store')
    parser.add_argument('-l', '--layer', dest='layer_idx', action='store',
                        type=int, help='layer index for plot')
    parser.add_argument('-p', '--plot', dest='plot3d', action='store_true',
                        help='plot the whole part')
    args = parser.parse_args()
    print(args)
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
    if args.plot3d:
        gcode_reader.plot()
    else:
        if args.layer_idx:
            gcode_reader.plot_layer(layer=args.layer_idx)
    # gcode_reader.plot_layers(min_layer=1, max_layer=2)
    # gcode_reader.plot()


if __name__ == "__main__":
    command_line_runner()

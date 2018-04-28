#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# Gcode reader
##################################

# standard library
import os.path
import sys
import argparse
import pprint

# third party library
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pp = pprint.PrettyPrinter(indent=4)


class GcodeReader:
    """ Gcode reader class """

    def __init__(self, filename, filetype="default"):
        if not os.path.exists(filename):
            print("{} does not exist!".format(filename))
            sys.exit(1)
        self.filename = filename
        self.filetype = filetype
        self.n_segs = 0
        self.segs = None
        self.n_layers = 0
        self.index_bars = []
        self._read()

    def _read(self):
        """ read the file and stores segs into data structure """
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

    def plot(self, color='blue'):
        """ plot the whole part in 3D """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        assert(self.n_segs > 0)
        x0, y0, x1, y1, z = self.segs[0, :]
        xs, ys, zs = [x0, x1], [y0, y1], [z, z]
        for x0, y0, x1, y1, z in self.segs[1:, :]:
            if x0 != xs[-1] or y0 != ys[-1] or z != zs[-1]:
                ax.plot(xs, ys, zs)
                xs, ys, zs = [x0, x1], [y0, y1], [z, z]
            else:
                xs.append(x1)
                ys.append(y1)
                zs.append(z)
        if len(xs) != 0:
            ax.plot(xs, ys, zs)
        plt.show()

    def plot_layer(self, layer=1):
        """ plot a specific layer in 2D """
        # make sure layer is in [1, self.n_layers]
        layer = max(layer, 1)
        layer = min(self.n_layers, layer)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        left, right = self.index_bars[layer - 1], self.index_bars[layer]
        x0, y0, x1, y1, _ = self.segs[left, :]
        xs, ys = [x0, x1], [y0, y1]
        for x0, y0, x1, y1, _ in self.segs[left + 1 : right, :]:
            if x0 != xs[-1] or y0 != ys[-1]:
                ax.plot(xs, ys)
                xs, ys = [x0, x1], [y0, y1]
            else:
                xs.append(x1)
                ys.append(y1)
        if len(xs) != 0:
            ax.plot(xs, ys)
        plt.show()


def command_line_runner():
    """ main function """
    parser = argparse.ArgumentParser(description='Gcode Reader')
    parser.add_argument(dest='gcode_file', help='gcode file', action='store')
    args = parser.parse_args()
    gcode_reader = GcodeReader(args.gcode_file)
    gcode_reader.plot_layer(layer=2)
    # gcode_reader.plot()


if __name__ == "__main__":
    command_line_runner()

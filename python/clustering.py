"""
Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella

This program is free software: you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as published
by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from read import read_assignments
import numpy as np
import os
import argparse


def contiguous_areas(ass):
    """ Returns areas covered by contiguous overlapping peptides from
    an assignment file """
    areas = []

    residues = [i + 1 for i in range(max(ass.T[2]))]
    coverage = np.zeros(max(ass.T[2]))
    for i in range(len(ass)):
        start = ass[i][1]
        end = ass[i][2] - 1
        for j in range(start, end):
            coverage[j] += 1

    limiting_residues = [residues[i] for i in np.where(coverage == 0)[0]]

    for i in range(len(limiting_residues) - 1):
        first = limiting_residues[i]
        second = limiting_residues[i + 1]
        length = second - first
        if length > 1:
            if i == len(limiting_residues) - 1:
                area = "%s-%s" % (first + 1, second + 1)
            else:
                area = "%s-%s" % (first + 1, second)
            areas.append(area)

    return areas


def run_mclust(areas):
    """ Calls the R script multi.r in folder R,
    which runs mclust clustering algortithm """
    for area in areas:
        txt = "Rscript ../R/multi.r "+area.split('-')[0]+' '+area.split('-')[1]
        os.system(txt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--ass")

    opts = parser.parse_args()

    if opts.ass:
        ass_file = opts.ass
        ass = read_assignments(ass_file)

    areas = contiguous_areas(ass)
    run_mclust(areas)

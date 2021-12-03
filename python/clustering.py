"""
Created on Wed 29 Sep 2021

@author: Michele Stofella
"""

from read import read_assignments
import numpy as np
import os
import argparse


def contiguous_areas(ass):
    """ Returns areas covered by contiguous overlapping peptides from
    an assignment file """
    areas = []

    residues = [i+1 for i in range(max(ass.T[2]))]
    coverage = np.zeros(max(ass.T[2]))
    for i in range(len(ass)):
        start = ass[i][1]
        end = ass[i][2]-1
        for j in range(start, end):
            coverage[j] += 1

    limiting_residues = [residues[i] for i in np.where(coverage == 0)[0]]

    for i in range(len(limiting_residues)-1):
        first = limiting_residues[i]
        second = limiting_residues[i+1]
        length = second-first
        if length > 1:
            if i == len(limiting_residues)-1:
                area = str(first+1)+'-'+str(second+1)
            else:
                area = str(first+1)+'-'+str(second)
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

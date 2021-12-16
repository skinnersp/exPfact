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

import os
import numpy as np
import pandas as pd
from read import read_pfact


def select_top_solutions(out_file, n):
    files = os.listdir()
    diff_files = []
    diff_value = []
    for file in files:
        if '.diff' in file and out_file in file:
            f = open(file, 'r')
            diff = pd.read_csv(file, header=None, delim_whitespace=True)
            diff_files.append(file)
            diff_value.append(np.average(diff[1]))
    sorted_files = pd.DataFrame(zip(diff_files, diff_value)).sort_values(1).reset_index(drop=True)
    np.savetxt('diff.list', sorted_files.values, fmt='%s %5.10f')

    p = []
    for i in range(int(sorted_files.shape[0]*n/100)):
        file = sorted_files[0][i]
        p.append(list(read_pfact(file.replace('.diff', '.pfact'))))

    with open("all.sp", "w+") as f:
        for i in range(len(p)):
            for j in range(len(p[i])):
                if j == len(p[i])-1:
                    f.write("{: <12}\n".format(round(p[i][j], 5)))
                else:
                    f.write("{: <12} ".format(round(p[i][j], 5)))


def run_descriptive():
    allsp = pd.read_csv("all.sp", header=None, delim_whitespace=True)

    with open('average.pfact', 'w+') as f:
        for i in range(allsp.shape[1]):
            round_mean = round(np.mean(allsp[i]), 5)
            round_std = round(np.std(allsp[i]), 5)
            f.write(str(i+1)+'\t'+str(round_mean)+'\t'+str(round_std)+'\n')

    with open('median.pfact', 'w+') as f:
        for i in range(allsp.shape[1]):
            round_median = round(np.median(allsp[i]), 5)
            f.write(str(i+1)+'\t'+str(round_median)+'\n')

    with open('minmax.pfact', 'w+') as f:
        for i in range(allsp.shape[1]):
            round_min = round(np.min(allsp[i]), 5)
            round_max = round(np.max(allsp[i]), 5)
            f.write(str(i+1)+'\t'+str(round_min)+'\t'+str(round_max)+'\n')


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--res")
    parser.add_argument("--top")

    opts = parser.parse_args()

    if opts.res:
        out_file = opts.res
    if opts.top:
        n = float(opts.top)
    else:
        n = 50

    select_top_solutions(out_file, n)
    run_descriptive()

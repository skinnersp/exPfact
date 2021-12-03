"""
Copyright (C) 2019-2020 Simon P. Skinner

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

import numpy as np


def write_pfact(params, fout_name):
    """
    Writes Pfactors to output file.
    :param params: array of pfactors
    :param fout_name: output file name.
    :return:
    """
    fout = open(fout_name + '.pfact', 'w')
    for ii, x in enumerate(params):
        fout.write("{} {}\n".format(ii + 1, x))
    fout.close()


def write_dpred(output_file, dpred, times, eps=0, suffix=".Dpred"):
    """
    Writes Dpred values to file.
    :param output_file: output file name
    :param dpred: array of Dpred values
    :param times: array of time points.
    :param suffix: suffix for output file.
    :return:
    """
    output_array = np.insert(dpred, [0], times, axis=0)
    if eps > 0:
        for i in range(1, len(output_array)):
            for j in range(len(output_array[i])):
                output_array[i][j] += np.random.normal(scale=eps)
                if output_array[i][j] > 1:
                    output_array[i][j] = 1
                elif output_array[i][j] < 0:
                    output_array[i][j] = 0
    np.savetxt(output_file + suffix, output_array.T, fmt='%.7g')


def write_diff(outfile, dpred, dexp):
    """
    Writes out root normalised differences between dpred and dexp
    :param outfile: output file name
    :param dpred: array of dpred values
    :param dexp: array of dexp values
    :return:
    """
    fout = open(outfile + '.diff', 'w')
    costs = [1/len(pred)*np.sum((pred-exp)**2) for pred, exp in zip(dpred, dexp)]
    for ii, cost in enumerate(costs):
        fout.write('{} {:e}\n'.format(ii + 1, cost))
    fout.close()


def write_combined_replicates(files, out):
    """
    Writes out .Dpred file where deuteration is mean of multiple replicates
    :param files: list of .Dpred files to be combined
    :param out:   name of output .Dpred file
    :return:
    """
    list_arrays = []
    weights = []
    for file in files:
        list_arrays.append(np.loadtxt(file))
        weights.append(np.loadtxt(file).T[1:])
    comb = np.mean(np.array(list_arrays), axis=0)

    all_w = []
    for i in range(len(weights[0])):
        for j in range(len(weights[0][i])):
            all_w.append([weights[k][i][j] for k in range(len(weights))])
    stds = [np.std(all_w[i]) for i in range(len(all_w))]
    pstd = np.sqrt(np.sum([std**2 for std in stds])/len(stds))

    print("  Pooled std: "+str(round(pstd, 5)))

    np.savetxt(out+'.Dpred', comb, fmt='%.7g')

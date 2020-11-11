"""
Copyright (C) 2019-2020 Simon P. Skinner

This program is free software: you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as published by
the Free Software Foundation.

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


def write_dpred(output_file, dpred, times, suffix=".Dpred"):
    """
    Writes Dpred values to file.
    :param output_file: output file name
    :param dpred: array of Dpred values
    :param times: array of time points.
    :param suffix: suffix for output file.
    :return:
    """
    output_array = np.insert(dpred, [0], times, axis=0)
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
    costs = [np.sqrt(1 / len(pred) * np.sum((pred - exp) ** 2)) for pred, exp in zip(dpred, dexp)]
    for ii, cost in enumerate(costs):
        fout.write('{} {}\n'.format(ii + 1, cost))
    fout.close()

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

########################################################################################################
#                                                                                                      #
#  This script fits a set of pfactors to a set of kint values from protein fragments.                  #
#                                                                                                      #
#  Compulsory arguments:                                                                               #
#                                                                                                      #
#  --temp: specifies the temperature for kint calculation.                                             #
#  --pH:   specifies the temperature for pH calculation.                                               #
#  --pfact: specifies file pfactors                                                                    #
#  --ass:  specifies file containing assignments of kints to dexp values.                              #
#  --times: list of time poinds at which D is calculated for each peptide in ass		               #
#  --seq: amino acid sequence file 								                                       #
#                                                                                                      #
#                                                                                                      #
#  Optional arguments:                                                                                 #
#  --base:  specifies directory where calculation files live (default: pwd)                            #
#  --out:   specifies name of output file for protection factors (default: pfact.out)                  #
#                                                                                                      #
########################################################################################################
"""

from calc_dpred import calculate_dpred
import os

from kint import calculate_kint_for_sequence

from read import read_assignments, \
        read_pfact, \
        read_seq,	 \
        read_time_points

from write import write_dpred, write_combined_replicates

import numpy as np

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--base")
parser.add_argument("--ass")
parser.add_argument("--pfact")
parser.add_argument("--times")
parser.add_argument("--temp")
parser.add_argument("--pH")
parser.add_argument("--seq")
parser.add_argument("--out")

parser.add_argument("--nrep")
parser.add_argument("--eps")

config = {}

opts = parser.parse_args()

# Compulsory arguments

if opts.base:
    config['base'] = opts.base
else:
    config['base'] = os.getcwd()

if opts.ass:
    config['assignments'] = opts.ass
if opts.temp:
    config['temperature'] = float(opts.temp)
if opts.pH:
    config['pH'] = float(opts.pH)
if opts.pfact:
    config['pfact'] = read_pfact(opts.pfact)
if opts.times:
    config['times'] = read_time_points(opts.times)
if opts.seq:
    config['sequence'] = read_seq(opts.seq)
    config['res1'] = 1
    config['resn'] = len(read_seq(opts.seq))

# Optional arguments
if opts.out:
    config['output'] = opts.out
else:
    config['output'] = None

if opts.nrep:
    nrep = int(opts.nrep)
else:
    nrep = 1
    
if opts.eps:
    eps = float(opts.eps)
else:
    eps = 0.
    
pfact = config['pfact']
assignments = read_assignments(config['assignments'])

assignment_set = set()
for ass in assignments:
    for x in range(int(ass[1]), int(ass[2]) + 1):
        assignment_set.add(x)

kint, prolines = calculate_kint_for_sequence(config['res1'], config['resn'], config['sequence'], config['temperature'], config['pH'])

outfiles = []
for i in range(1,nrep+1):
    if nrep > 1:
        outfile = config['output']+str(i)
    else:
        outfile = config['output']

    dpred = calculate_dpred(pfact, config['times'], kint, assignments)
    write_dpred(outfile, dpred, config['times'], eps)
    outfiles.append(outfile+".Dpred")
    
if nrep > 1:
    write_combined_replicates(outfiles, out=config['output'])

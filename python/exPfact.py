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

###############################################################################
#
#  This script fits a set of pfactors to a set of kint values
#  from protein fragments.
#
#  Compulsory arguments:
#
#  --temp: specifies the temperature for kint calculation.
#  --pH:   specifies the temperature for pH calculation.
#  --dexp: specifies file containiing dexp values.
#  --ass:  specifies file containing assignments of kints to dexp values.
#
#
#  Optional arguments:
#  --base:  specifies directory where calculation files live (default: pwd)
#  --out:   specifies name of output file for protection factors
#           (default: pfact.out)
#  --rand:  specifies number of iterations of montecarlo steps to be performed
#           to determine initial pfactors (default: 1000)
#  --tol:   tolerance threshold for least squares convergence (default: 1e-6)
#  --harm:  introduces a penalty for large differences in predicted lnP
#           for adjacent residues
#
###############################################################################
"""

import os
import sys
import argparse

from calc_dpred import calculate_dpred

from calculate import cost_function,\
    do_random_search, \
    fit_pfact

from kint import calculate_kint_for_sequence

from read import read_assignments, \
    read_configuration, \
    read_dexp, \
    read_pfact, \
    read_seq

from write import write_diff, \
    write_dpred, \
    write_pfact

from logger import log

def run(base_dir, dexp, assignments, pfact, random_steps, time_points,
        harmonic_term, output_file, tolerance, weights, pH, temperature,
        seq, res1, resn):
    """

    :param base_dir: base directory for all input files.
    :param dexp: file containing dexp values.
    :param assignments: file containing assignments of kints to dexp values.
    :param pfact: file containing pfactor values.
    :param random_steps: number of steps for random search.
    :param time_points: a list of experiment time points.
    :param harmonic_term: term to be used for harmonic cost scoring.
    :param output_file: stub for all output files.
    :param tolerance: tolerance value for minimisation convergence.
    :return:
    """

    assignment_set = set()
    for ass in assignments:
        for x in range(int(ass[1]), int(ass[2]) + 1):
            assignment_set.add(x)

    pfactor_filter = set()
    for ass in assignments:
        for x in range(int(ass[1] + 1), int(ass[2]) + 1):
            pfactor_filter.add(x)
        if ass[1] < min(pfactor_filter):
            pfactor_filter.add(ass[1])

    kint, prolines = calculate_kint_for_sequence(res1, resn, seq,
                                                 temperature, pH)

    if not pfact:
        if random_steps:
            rand_output = do_random_search(
                                           kint,
                                           random_steps,
                                           pfactor_filter,
                                           dexp,
                                           time_points,
                                           assignments,
                                           harmonic_term,
                                           prolines,
                                           weights,
                                           seed=None
                                           )
            min_score = min(rand_output.keys())
            init_array = rand_output[min_score]
        else:
            init_array = [1 if ii not in prolines or ii == 0
                          or ii + 1 in pfactor_filter
                          else -1 for ii in range(max(pfactor_filter))]

    else:
        init_array = read_pfact(pfact)

    bounds = [(0.00001, 20) if x >= 0 else (-1, -1) if x == -1
              else (0, 0) for x in init_array]

    pfit = fit_pfact(init_array, dexp, time_points, assignments,
                     harmonic_term, kint, bounds, tolerance, weights)

    write_pfact(pfit.x, output_file)

    dpred = calculate_dpred(pfit.x, time_points, kint, assignments)

    write_dpred(output_file, dpred, time_points)
    write_diff(output_file, dpred, dexp)

    final_score = cost_function(pfit.x, dexp, time_points, assignments,
                                harmonic_term, kint, weights)
    print('Final value of cost function w harm term: {}'.format(final_score))
    final_score = cost_function(pfit.x, dexp, time_points,
                                assignments, 0.0, kint, weights)
    print('Final value of cost function w/o harm term: {}'.format(final_score))


def main(argv):
    """
    :param argv: input arguments from command line.
    :return:
    """
    log.info("Running exPfact.py")

    parser = argparse.ArgumentParser()

    parser.add_argument("--base")
    parser.add_argument("--dexp")
    parser.add_argument("--ass")
    parser.add_argument("--pfact")
    parser.add_argument("--weights")
    parser.add_argument("--out")
    parser.add_argument("--predict")
    parser.add_argument("--times")
    parser.add_argument("--tol")
    parser.add_argument("--harm")
    parser.add_argument("--rand")
    parser.add_argument("--temp")
    parser.add_argument("--pH")
    parser.add_argument("--seq")
    parser.add_argument("--rep")

    if sys.argv[1].endswith('.json'):
        config = read_configuration(sys.argv[1])
    else:

        config = {}

        opts = parser.parse_args()

        # Compulsory arguments

        if opts.base:
            config['base'] = opts.base
            print("Base directory= ", config['base'])
        else:
            config['base'] = os.getcwd()

        if opts.dexp:
            config['dexp'], config['times'] = read_dexp(opts.dexp)
        if opts.ass:
            config['assignments'] = opts.ass
            print("ass= ", config['assignments'])
        if opts.temp:
            config['temperature'] = float(opts.temp)
        if opts.pH:
            config['pH'] = float(opts.pH)

        if opts.seq:
            config['sequence'] = read_seq(opts.seq)
            config['res1'] = 1
            config['resn'] = len(read_seq(opts.seq))

        # Optional arguments

        if opts.predict:
            config['predict'] = True

        if opts.times:
            config['time_points'] = opts.times

        if opts.pfact:
            config['pfact'] = opts.pfact
            print("pfactfile= ", config['pfact'])
        else:
            config['pfact'] = None

        if opts.out:
            config['output'] = opts.out
        else:
            config['output'] = None

        if opts.rand:
            config['do_random_search'] = True
            config['random_search_steps'] = int(opts.rand)
        else:
            config['do_random_search'] = False
            config['random_search_steps'] = None
        if opts.tol:
            config['tolerance'] = float(opts.tol)
        else:
            config['tolerance'] = None
        if opts.harm:
            config['harmonic_factor'] = float(opts.harm)
        else:
            config['harmonic_factor'] = 0

        if opts.weights:
            config['weights'] = read_dexp(opts.weights)[0]
        else:
            config['weights'] = None

        if opts.rep:
            n_rep = int(opts.rep)
        else:
            n_rep = 1

    assignments = read_assignments(config['assignments'])

    for i in range(n_rep):
        log.info("Minimization %s of %s" % (str(i),str(n_rep)))
        if n_rep > 1:
            outfile = config['output']+str(i)
        else:
            outfile = config['output']

        run(
            config['base'],
            config['dexp'],
            assignments,
            config['pfact'],
            config['random_search_steps'],
            config['times'],
            config['harmonic_factor'],
            outfile,
            config['tolerance'],
            config['weights'],
            config['pH'],
            config['temperature'],
            config['sequence'],
            config['res1'],
            config['resn']
            )


if __name__ == "__main__":
    try:
        sys.argv[1]
    except IndexError:
        print(__doc__)
        log.error("Error while running ExPfact.py")
        exit()
    main(sys.argv[1:])

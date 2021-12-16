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

import pandas as pd
import numpy as np
import glob
import re
from sklearn.metrics import r2_score
from itertools import combinations

from read import read_assignments, read_seq, read_pfact, read_time_points
from kint import calculate_kint_for_sequence
from kback import calculate_kback_for_sequence
from Hisotope import fully_protonated_envelope


def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split('([0-9]+)', s)]


def natural_sort(lst):
    """ Sort the given list in the way that humans expect.
    """
    lst.sort(key=alphanum_key)


def deuterium_uptake(t, kint, lnp):
    '''Calculates deuterium uptake at time t given kint and lnP'''
    p = np.exp(lnp)
    tmp = 0
    namide = 0
    for i in range(0, len(kint)):
        if kint[i] >= 0:
            namide += 1
            tmp += np.exp(-kint[i] / p[i] * t * 60)
    return (namide - tmp)/namide


def single_deuterium_uptake(t, kint, lnp):
    '''Calculates deuterium uptake at time t given kint and lnP'''
    p = np.exp(lnp)
    if kint > 0:
        return 1 - np.exp(-kint / p * t * 60)
    else:
        return 0


def single_deuterium_downtake(t, kback, lnp):
    '''Calculates deuterium downtake at time t given kback and lnP'''
    p = np.exp(lnp)
    if kback > 0:
        return 1 - np.exp(-kback / p * t * 60)
    else:
        return 0


def set_a_group(a, k):
    return list(combinations(a, k))


def probability_uptake(k, t, kint, lnp):
    ''' Calculates the probability that k residues have exchanged at time t
    given intrinsic rates kint and protection factors lnP '''
    a = set_a_group([i for i in range(len(kint))], k)
    probability = 0
    for element in a:
        d1 = []
        d2 = []
        for j in range(len(kint)):
            if j in element:
                d1.append(single_deuterium_uptake(t, kint[j], lnp[j]))
            else:
                d2.append(1 - single_deuterium_uptake(t, kint[j], lnp[j]))
        p1 = np.prod(d1)
        p2 = np.prod(d2)
        probability += p1 * p2
    return probability


def probability_downtake(k, t, kint, lnp):
    ''' Calculates the probability that k residues have exchanged at time t
    given intrinsic rates kint and protection factors lnP '''
    a = set_a_group([i for i in range(len(kint))], k)
    probability = 0
    for element in a:
        d1 = []
        d2 = []
        for j in range(len(kint)):
            if j in element:
                d1.append(single_deuterium_downtake(t, kint[j], lnp[j]))
            else:
                d2.append(1-single_deuterium_downtake(t, kint[j], lnp[j]))
        p1 = np.prod(d1)
        p2 = np.prod(d2)
        probability += p1 * p2
    return probability


def isotopic_envelope(t, kint, lnp, exchange):
    """ Returns the probability intensities of isotopic envelope at time t for

    * exchange == 'f': forward exchange (protonated protein in D2O)
    * exchange == 'b': back exchange(deuterated protein in H2O)
    """
    k_values = [i for i in range(2*len(kint))]
    freq = []
    for k in k_values:
        if exchange == 'f':
            freq.append(probability_uptake(k, t, kint, lnp))
        elif exchange == 'b':
            freq.append(probability_downtake(k, t, kint, lnp))
    return freq


def centered_isotopic_envelope(t, kint, lnp, fr0):
    fr = isotopic_envelope(t, kint, lnp, exchange='f')
    f = np.zeros(2 * len(kint))
    for i in range(len(f)):
        for j in range(i + 1):
            f[i] += fr0[i-j] * fr[j]
    return f


def back_centered_isotopic_envelope(t, kint, lnp, fr0):
    fr = isotopic_envelope(t, kint, lnp, exchange='b')
    f = np.zeros(2 * len(kint))
    for i in range(len(f) + 1, -1, -1):
        for j in range(len(f) - i):
            f[i] += fr0[i + j] * fr[j]
    f = [f[i] / sum(f) * 100 for i in range(len(f))]
    return f


def predict_isotopic_envelope(ass_file, seq_file, temperature, pH,
                              lnp_file, times_file, pep, charge_state,
                              exchange, out_file, pi0_file=''):

    seq = read_seq(seq_file)
    times = read_time_points(times_file)

    # Select residues involving the selected peptide
    ass = read_assignments(ass_file)
    start_res = ass[int(pep) - 1][1]
    end_res = ass[int(pep) - 1][2]

    # Upload kint and lnP values
    if exchange == 'f':
        kint, _ = calculate_kint_for_sequence(1, len(seq), seq,
                                              float(temperature), float(pH))
        kint = kint[start_res:end_res]
    elif exchange == 'b':
        kint, _ = calculate_kback_for_sequence(1, len(seq), seq,
                                               float(temperature), float(pH))
        kint = kint[start_res:end_res]

    lnP = read_pfact(lnp_file)[start_res:end_res]
    # Calculate fully protonated isotopic envelope
    if exchange == 'f':
        pi0 = fully_protonated_envelope(seq[start_res:end_res + 1],
                                        z=charge_state)
        mass = list(pi0.keys())
        fr0 = list(pi0.values())
        while len(mass) <= 2 * len(kint[start_res:end_res + 1]):
            mass.append((mass[-1] + 1.00627 * int(charge_state))/charge_state)
            fr0.append(0)
            print(mass, fr0)
    elif exchange == 'b':
        pi0 = pd.read_csv(pi0_file, skiprows=1,
                          header=None, delim_whitespace=True)
        mass = list(pi0[1])
        u_fr0 = list(pi0[2])
        fr0 = centered_isotopic_envelope(0, kint, lnP, u_fr0)

    # Calculate isotopic envelopes at different times
    for i in range(len(times)):
        if exchange == 'f':
            f1 = centered_isotopic_envelope(times[i], kint, lnP, fr0)
        elif exchange == 'b':
            f1 = back_centered_isotopic_envelope(times[i], kint, lnP, fr0)

        f1 = [f1[j] / sum(f1) * 100 for j in range(len(f1))]
        with open("%s.%s.isot" % (out_file, str(i)), 'w+') as f:
            f.write('# '+seq[start_res:end_res]+'\n')
            for j in range(len(f1)):
                f.write('%d\t' % j)
                f.write('%5.5f\t' % mass[j])
                f.write('%5.2f\t' % f1[j])
                last_col = f1[j] / max(f1) * 100
                if j == len(f1) - 1:
                    f.write('%5.2f' % last_col)
                else:
                    f.write('%5.2f\n' % last_col)


def generate_back_exchange_time_points(start=-5, end=5, num=500):
    with open('back.times', 'w+') as f:
        a = np.logspace(start, end, num)
        for i in range(0, len(a)):
            f.write('%5.15f\n' % a[i])


def sticks_from_experimental_envelope(exp_env, corr_env, z):

    mass = list(corr_env[1])
    fr = np.zeros(len(mass))
    for i in range(len(exp_env)):
        for j in range(len(mass)):
            if exp_env[0][i] * z >= mass[j] - 0.5 and exp_env[0][i] * z < mass[j] + 0.5:
                fr[j] += exp_env[1][i]
    fr = fr / sum(fr) * 100

    return mass, fr


def compare_predictions(fr, prefix):

    times = read_time_points('back.times')

    files = []
    for file in glob.glob(prefix):
        files.append(file)
    sort_nicely(files)

    scores = []
    for i in range(len(files)):

        back_env = pd.read_csv(files[i], header=None, sep='\t', skiprows=1)

        score = r2_score(fr, back_env[2])
        scores.append(score)

    res = pd.DataFrame([files, times, scores]).transpose()

    return res

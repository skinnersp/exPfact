"""
Created on Thu 02 Dec 2021
Unit tests

@author: Michele Stofella
"""

import pytest
import numpy as np
import math

from kint import calculate_kint_for_sequence
from kback import calculate_kback_for_sequence
from calculate import calculate_rms, \
                      do_random_search, \
                      predict_dexp
from Hisotope import fully_protonated_envelope
from isenv_functions import isotopic_envelope

##############################################################################
#
# Testing intrinsic exchange rates calculations
# Scripts: kint.py
#          kback.py
#
##############################################################################

@pytest.mark.parametrize("seq, kint_englander",
    [("SAMPLE", np.array([-1.0, 1.3e+06, 1.6e+04, -1.0, 2.4e+03, 1.2e+02])),
     ("SICILY", np.array([-1.0, 2.5e+05, 3.9e+04, 1.1e+04, 2.5e+03, 8.4e+01])),
     ("RIVENDELL", np.array([-1.0, 2.1e+05, 1.9e+03, 9.1e+03, 3.6e+04,
                             4.2e+04, 8.4e+03, 3.1e+03, 4.1e+01]))
    ])
def test_forward_intrinsic_rates(seq, kint_englander):
    """ Checks that forward intrinsic exchange rates are correctly calculated
    by the script kint.py. The results are tested against the rates obtained
    for the same sequence by the Englander group excel spreadsheet """
    # intrinsic exchange rates calculated by kint.py
    kint, pro = calculate_kint_for_sequence(1,len(seq),seq,300,7)
    for i in range(len(seq)):
        # check that the rates are the same (maximum difference 1%)
        assert np.abs(kint[i]/kint_englander[i]-1) < 1


@pytest.mark.parametrize("seq, kback_englander",
    [("SAMPLE", np.array([-1.0, 6.7e+06, 7.9e+04, -1.0, 1.2e+04, 6.1e+02])),
     ("SICILY", np.array([-1.0, 1.2e+06, 2.0e+06, 5.3e+04, 1.3e+04, 4.2e+02])),
     ("RIVENDELL", np.array([-1.0, 1.0e+06, 9.4e+03, 4.5e+04, 1.8e+05,
                             2.1e+05, 4.1e+04, 1.5e+04, 2.1e+02]))
    ])
def test_backward_intrinsic_rates(seq, kback_englander):
    """ Checks that backward intrinsic exchange rates are correctly calculated
    by the script kback.py. The results are tested against the rates obtained
    for the same sequence by the Englander group excel spreadsheet """
    # intrinsic exchange rates calculated by kback.py
    kback, pro = calculate_kback_for_sequence(1,len(seq),seq,300,7)
    for i in range(len(seq)):
        assert np.abs(kback[i]/kback_englander[i]-1) < 1


@pytest.mark.parametrize("seq, expected_prolines",
    [("AAPAAPAAPAA", [3, 6, 9]),
     ("APAPAPAPAPA", [2, 4, 6, 8, 10])])
def test_prolines_kint(seq, expected_prolines):
    """ Check that the script kint.py correctly identifies prolines along
    the sequence of the peptide and that the intrinsic exchange rate at those
    residue is set to -1.0 """
    kint, prolines = calculate_kint_for_sequence(1,len(seq),seq,300,7)
    assert len(prolines) == len(expected_prolines)
    for i in range(len(prolines)):
        assert prolines[i] == expected_prolines[i]
        assert kint[prolines[i]-1] < 0


@pytest.mark.parametrize("seq, expected_prolines",
    [("AAPAAPAAPAA", [3, 6, 9]),
     ("APAPAPAPAPA", [2, 4, 6, 8, 10])])
def test_prolines_kback(seq, expected_prolines):
    """ Check that the script kback.py correctly identifies prolines along
    the sequence of the peptide and that the intrinsic exchange rate at those
    residue is set to -1.0 """
    kback, prolines = calculate_kback_for_sequence(1,len(seq),seq,300,7)
    assert len(prolines) == len(expected_prolines)
    for i in range(len(prolines)):
        assert prolines[i] == expected_prolines[i]
        assert kback[prolines[i]-1] < 0


##############################################################################
#
# Testing cost function calculations
# Scripts: calculate.py
#
##############################################################################


@pytest.mark.parametrize("dpred, dexp, nj, weights, rms",
    [(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
      np.array([[1.1, 2.1, 3.1], [4.2, 4.9, 5.8], [7.3, 8.9, 9.4]]),
      12,
      None,
      0.09833),
     (np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
      np.array([[1.1, 2.1, 3.1], [4.2, 4.9, 5.8], [7.3, 8.9, 9.4]]),
      15,
      np.array([[1, 2, 1], [1, 1, 1], [0.5, 0.7, 1]]),
      0.06013),
     (np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
      np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
      5,
      None,
      0.00000)
      ])
def test_calculate_rms(dpred, dexp, nj, weights, rms):
    """ Checks that the root mean square calculation is correct for cases
    calculated manually and that it is zero when a dataset is compared to
    itself """
    r = calculate_rms(dpred, dexp, nj, weights)
    assert np.round(r, 5) == rms


@pytest.mark.parametrize("kint, search_steps, pfactor_filter, dexp, time_points, assignments, harmonic_term, prolines, weights, seed",
    [(np.array([-1.0, 5.2e+02, 3.2e+04, 3.5e+04, 9.2e+05,
                -1.0, 7.2e+03, 2.2e+02, 8.1e+02, 4.0e+01]),
     10,
     {1,2,3,4,5,6,7,8,9,10},
     np.array([[.10, .25, .48], [.31, .48, .77], [.21, .58, .88],
               [.15, .50, .91], [.58, .79, .95]]),
     np.array([0.00139, 0.08333, 0.7500]),
     np.array([[1, 1, 6], [2, 1, 10], [3, 3, 8], [4, 5, 9], [5, 5, 10]]),
     0,
     [6],
     None,
     31415)])
def test_random_search(kint, search_steps, pfactor_filter, dexp, time_points,
                       assignments, harmonic_term, prolines, weights, seed):
    """ Checks that the random search is repeated the set number of times,
    that prolines are correclty identified and set to -1.0 """
    score_array = do_random_search(kint, search_steps, pfactor_filter, dexp,
                                   time_points, assignments,
                                   harmonic_term, prolines, weights, seed)

    assert len(score_array) == search_steps
    for i in range(len(score_array)):
        assert score_array[list(score_array.keys())[i]][5] < 0


@pytest.mark.parametrize("pfact, time_points, kint, assignments, expected", [
    (np.array([-1., 3., 2., 5., 1.]),
     np.array([0., 99999.]),
     np.array([-1., 5e+03, 2e+05, 5e+03, 1e+02]),
     np.array([[1, 1, 5], [2, 3, 5]]),
     np.array([[0., 1.], [0., 1.]]))
])
def test_predict_dexp(pfact, time_points, kint, assignments, expected):
    """ Checks that the deuterium uptake is 0 at time 0 and 1 at long times """
    dpred = predict_dexp(pfact, time_points, kint, assignments)
    np.testing.assert_array_equal(dpred, expected)


##############################################################################
#
# Testing isotopic envelope calculations
# Scripts: Hisotope.py
#          isenv.py
#          isenv_functions.py
#
##############################################################################


@pytest.mark.parametrize("sequence, z, expected", [
    ("SAMPLE",
    3,
    {216.4: 66.71,
     216.8: 22.88,
     217.1: 8.12,
     217.4: 1.87,
     217.8: 0.36}),
    ("SICILY",
    1,
    {711.4: 62.47,
     712.4: 25.68,
     713.4: 9.06,
     714.4: 2.26,
     715.4: 0.45,
     716.4: 0.07,
     717.4: 0.01}),
    ("RIVENDELL",
    2,
    {550.8: 53.67,
     551.3: 31.64,
     551.8: 11.10,
     552.3: 2.86,
     552.8: 0.60,
     553.3: 0.11,
     553.8: 0.02})
    ])
def test_fully_protonated_envelope(sequence, z, expected):
    """ Checks the accuracy of the calculation of the fully protonated envelope;
    the results are parametrized following the online software MS-Isotope
    (link: https://prospector.ucsf.edu/prospector/) """
    env = fully_protonated_envelope(sequence, z, write=False)
    for i in range(len(expected)):
        assert list(expected.keys())[i] == round(list(env.keys())[i], 1)
        assert np.abs(list(expected.values())[i]-list(env.values())[i]) < 1


def test_isotopic_envelope(t, kint, lnP, exchange):

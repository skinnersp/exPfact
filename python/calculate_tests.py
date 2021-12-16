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

import pytest
import numpy as np

from calculate import calculate_rms, \
                      do_random_search, \
                      predict_dexp


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

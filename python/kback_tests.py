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
import math

from kback import calculate_kback_for_sequence

##############################################################################
#
# Testing intrinsic exchange rates calculations
# Scripts: kback.py
#
##############################################################################

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
def test_prolines_kback(seq, expected_prolines):
    """ Check that the script kback.py correctly identifies prolines along
    the sequence of the peptide and that the intrinsic exchange rate at those
    residue is set to -1.0 """
    kback, prolines = calculate_kback_for_sequence(1,len(seq),seq,300,7)
    assert len(prolines) == len(expected_prolines)
    for i in range(len(prolines)):
        assert prolines[i] == expected_prolines[i]
        assert kback[prolines[i]-1] < 0

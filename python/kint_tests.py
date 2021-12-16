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

from kint import calculate_kint_for_sequence

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

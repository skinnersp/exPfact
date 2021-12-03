"""
Created on Thu 02 Dec 2021
Unit tests

@author: Michele Stofella
"""

import pytest
from hypothesis import given
import hypothesis.strategies as st

from kint import calculate_kint_for_sequence
from kback import calculate_kback_for_sequence
import numpy as np
import math

def test_forward_intrinsic_rates():
    """ Checks that forward intrinsic exchange rates are correctly calculated
    by the script kint.py. The results are tested against the rates obtained
    for the same sequence by the Englander group excel spreadsheet """
    seq = 'SAMPLE'
    # intrinsic exchange rates calculated by kint.py
    kint, pro = calculate_kint_for_sequence(1,len(seq),seq,300,7)
    # intrinsic exchange rates calculated by Englander group excel spreadsheet
    kint_englander = [-1.0, 1.3e+06, 1.6e+04, -1.0, 2.4e+03, 1.2e+02]
    for i in range(len(seq)):
        # check that the rates are the same (maximum difference 1%)
        assert np.abs(kint[i]/kint_englander[i]-1) < 1


def test_backward_intrinsic_rates():
    """ Checks that backward intrinsic exchange rates are correctly calculated
    by the script kback.py. The results are tested against the rates obtained
    for the same sequence by the Englander group excel spreadsheet """
    seq = 'SAMPLE'
    # intrinsic exchange rates calculated by kback.py
    kback, pro = calculate_kback_for_sequence(1,len(seq),seq,300,7)
    # intrinsic exchange rates calculated by Englander group excel spreadsheet
    kback_englander = [-1.0, 6.7e+06, 7.9e+04, -1.0, 1.2e+04, 6.1e+02]
    for i in range(len(seq)):
        assert np.abs(kback[i]/kback_englander[i]-1) < 1


def test_prolines_kint():
    """ Check that the script kint.py correctly identifies prolines along
    the sequence of the peptide and that the intrinsic exchange rate at those
    residue is set to -1.0 """
    seq = 'AAPAAPAAPAA'
    kint, prolines = calculate_kint_for_sequence(1,len(seq),seq,300,7)
    assert len(prolines) == 3
    for i in range(len(prolines)):
        assert kint[prolines[i]-1] < 0


def test_prolines_kback():
    """ Check that the script kback.py correctly identifies prolines along
    the sequence of the peptide and that the intrinsic exchange rate at those
    residue is set to -1.0 """
    seq = 'AAPAAPAAPAA'
    kback, prolines = calculate_kback_for_sequence(1,len(seq),seq,300,7)
    assert len(prolines) == 3
    for i in range(len(prolines)):
        assert kback[prolines[i]-1] < 0


if __name__ == '__main__':
    seq = 'AAPAAPAAPAA'
    kint, prolines = calculate_kint_for_sequence(1,len(seq),seq,300,7)
    print(prolines)
    print(kint)

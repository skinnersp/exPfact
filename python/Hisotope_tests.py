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

from Hisotope import fully_protonated_envelope

##############################################################################
#
# Testing isotopic envelope calculations
# Scripts: Hisotope.py
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

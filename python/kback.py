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

from math import log10
import constants_DH as cst
import numpy as np


def calculate_kback_for_sequence(first_residue, last_residue, seq, temperature, pH):
    """
    Calculates the intrinsic exchange rates for a deuterated protein/peptide in water.
    Based on the work of Bai et al. (1993), Connelly et al. (1993), Nguyen et al. (2018)

    Parameters
    ----------
    first_residue : int
        Residue index along the sequence of the protein of the first amino
        acid of the sequence given as input.
    last_residue : int
        Residue index along the sequence of the protein of the last amino
        acid of the sequence given as input.
    seq : str
        Amino acid sequence (letter code) of the protein/peptide.
    temperature : float
        Temperature of the solution (in kelvin).
    pH : float
        Corrected pD: pDcorr = pDread+0.4

    Returns
    -------
    kback : ndarray
        Intrinsic exchange rates of the deuterated peptide in H2O in hr-1.
        The rate of prolines is set to -1.0.
    prolines : list
        List containing residue index of prolines.

    """
    prolines = []
    kback = np.zeros((last_residue))
    kback.fill(-1)
    res1 = ""
    jj = 0
    for assignment in range(1, len(seq)+1):
        res = seq[jj]
        if not res1 == "":
            if assignment - first_residue == 0:
                kback[assignment - 1] = -1
            elif seq[jj] == "P":
                kback[assignment - 1] = -1
                prolines.append(first_residue + jj)
            else:
                kback[assignment - 1] = calculate_kback_per_residue(res1, res,
                                                                  assignment,
                                                                  len(seq),
                                                                  temperature,
                                                                  pH)
        jj += 1
        res1 = res

    return kback, prolines


def calculate_kback_per_residue(residue1, residue2, num, length,
                                temperature, pH):
    """
    This function calculates the kback of a residue.
    The first argument is the residue i and the second the residue i-1.
    """

    lamb1 = acid(residue2, temperature, pH, "lamb")
    rho1 = acid(residue1, temperature, pH, "rho")
    if num == 2:
        rho1 += cst.rho_Nterm_acid
    elif (num == length):
        lamb1 += cst.lamb_Cterm_acid
    Fa = 10**(lamb1+rho1)

    lamb2 = base(residue2, temperature, pH, "lamb")
    rho2 = base(residue1, temperature, pH, "rho")
    if num == 2:
        rho2 += cst.rho_Nterm_base
    elif (num == length):
        lamb2 += cst.lamb_Cterm_base
    Fb = 10**(lamb2+rho2)

    kback = Fa * cst.ka * cst.get_D(pH) * cst.get_Fta(temperature) * 3600 +\
        Fb * cst.kb * cst.get_OD(pH) * cst.get_Ftb(temperature) * 3600 +\
        Fb * cst.kw * cst.get_Ftw(temperature) * 3600

    return kback


def acid(residue, temperature, pH, value):

    if residue == "H":
        lamb = log10(10**(-0.80-pH)/(10**(-cst.get_pK_his(temperature))+10**(-pH))+10**(0.00-cst.get_pK_his(temperature))/(10**(-cst.get_pK_his(temperature))+10**(-pH)))
        rho = log10(10**(-0.51-pH)/(10**(-cst.get_pK_his(temperature))+10**(-pH))+10**(0.00-cst.get_pK_his(temperature))/(10**(-cst.get_pK_his(temperature))+10**(-pH)))
    elif residue == "D":
        lamb = log10(10**(-0.90-pH)/(10**(-cst.get_pK_asp(temperature))+10**(-pH))+10**(0.90-cst.get_pK_asp(temperature))/(10**(-cst.get_pK_asp(temperature))+10**(-pH)))
        rho = log10(10**(-0.12-pH)/(10**(-cst.get_pK_asp(temperature))+10**(-pH))+10**(0.58-cst.get_pK_asp(temperature))/(10**(-cst.get_pK_asp(temperature))+10**(-pH)))
    elif residue == "E":
        lamb = log10(10**(-0.60-pH)/(10**(-cst.get_pK_glu(temperature))+10**(-pH))+10**(-0.90-cst.get_pK_glu(temperature))/(10**(-cst.get_pK_glu(temperature))+10**(-pH)))
        rho = log10(10**(-0.27-pH)/(10**(-cst.get_pK_glu(temperature))+10**(-pH))+10**(0.31-cst.get_pK_glu(temperature))/(10**(-cst.get_pK_glu(temperature))+10**(-pH)))
    else:
        lamb = cst.para[residue][0]
        rho = cst.para[residue][1]
    if value == "lamb":
        return lamb
    elif value == "rho":
        return rho


def base(residue, temperature, pH, value):

    if residue == "H":
        lamb = log10(10**(0.80-pH)/(10**(-cst.get_pK_his(temperature))+10**(-pH))+10**(-0.10-cst.get_pK_his(temperature))/(10**(-cst.get_pK_his(temperature))+10**(-pH)))
        rho = log10(10**(0.83-pH)/(10**(-cst.get_pK_his(temperature))+10**(-pH))+10**(0.14-cst.get_pK_his(temperature))/(10**(-cst.get_pK_his(temperature))+10**(-pH)))
    elif residue == "D":
        lamb = log10(10**(0.69-pH)/(10**(-cst.get_pK_asp(temperature))+10**(-pH))+10**(0.10-cst.get_pK_asp(temperature))/(10**(-cst.get_pK_asp(temperature))+10**(-pH)))
        rho = log10(10**(0.60-pH)/(10**(-cst.get_pK_asp(temperature))+10**(-pH))+10**(-0.18-cst.get_pK_asp(temperature))/(10**(-cst.get_pK_asp(temperature))+10**(-pH)))
    elif residue == "E":
        lamb = log10(10**(0.24-pH)/(10**(-cst.get_pK_glu(temperature))+10**(-pH))+10**(-0.11-cst.get_pK_glu(temperature))/(10**(-cst.get_pK_glu(temperature))+10**(-pH)))
        rho = log10(10**(0.39-pH)/(10**(-cst.get_pK_glu(temperature))+10**(-pH))+10**(-0.15-cst.get_pK_glu(temperature))/(10**(-cst.get_pK_glu(temperature))+10**(-pH)))
    else:
        lamb = cst.para[residue][2]
        rho = cst.para[residue][3]

    if value == "lamb":
        return lamb
    elif value == "rho":
        return rho

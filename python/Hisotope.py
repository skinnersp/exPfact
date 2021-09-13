"""
Created on Fri 13 Sep 2021

@author: Michele Stofella
"""

from pyopenms import *

def fully_protonated_envelope(sequence, z):
    """
    Calculates the fully protonated envelope of a peptide with given
    sequence and charge state. 

    Parameters
    ----------
    sequence : str
        Sequence of the peptide (letter code)
    z : int
        Charge state of the peptide.

    Returns
    -------
    isenv : dict
        Dictionary encoding the m/z values of the isotopic envelope. 
        The length of the dictionary is twice the length of the sequence. 
    
    The script also generates a 'sequence.txt' file containing the 
    predicted isotopic envelope in the format:
        isotope_index m/z intensity(sum=100) intensity(max=100)
    """
    seq = AASequence.fromString(sequence)
    seq_formula = seq.getFormula()

    isenv = {}
    isotopes = seq_formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(2*len(sequence)))
    for iso in isotopes.getContainer():
        isenv[(iso.getMZ()+z)/z] = iso.getIntensity()*100
    
    with open(sequence+'.txt','w+') as f:
        f.write("# "+sequence+"\n")
        for i in range(len(isenv)):
            isotope = list(isenv.keys())[i]
            intens1 = isenv[isotope]
            intens2 = isenv[isotope]/max(list(isenv.values()))*100
            f.write("%d %5.5f %5.2f %5.2f\n" % (i, isotope, intens1, intens2))
    print("Fully protonated envelope saved in file "+sequence+".txt!")
    
    return isenv
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

from pyopenms import AASequence, \
                     CoarseIsotopePatternGenerator


def fully_protonated_envelope(sequence, z, write=True):
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
    isotopic_envelope : dict
        Dictionary encoding the m/z values of the isotopic envelope.
        The length of the dictionary is twice the length of the sequence.

    The script also generates a 'sequence.txt' file containing the
    predicted isotopic envelope in the format:
        isotope_index m/z intensity(sum=100) intensity(max=100)
    """
    seq = AASequence.fromString(sequence)
    seq_formula = seq.getFormula()

    isotopic_envelope = {}
    isotopes = seq_formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(2*len(sequence)))
    for iso in isotopes.getContainer():
        isotopic_envelope[(iso.getMZ()+z)/z] = iso.getIntensity()*100

    if write:
        with open(sequence+'.txt', 'w+') as f:
            f.write("# "+sequence+"\n")
            for i in range(len(isotopic_envelope)):
                isotope = list(isotopic_envelope.keys())[i]
                intens1 = isotopic_envelope[isotope]
                intens2 = isotopic_envelope[isotope]/max(list(isotopic_envelope.values()))*100
                f.write("%d %5.5f %5.2f %5.2f\n" % (i, isotope, intens1, intens2))
        print("Fully protonated envelope saved in file "+sequence+".txt!")

    return isotopic_envelope


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--seq")
    parser.add_argument("--z")

    opts = parser.parse_args()

    if opts.seq:
        sequence = opts.seq

    if opts.z:
        z = int(opts.z)

    fully_protonated_envelope(sequence, z)

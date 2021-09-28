from pyteomics.parser import cleave, expasy_rules
from pyteomics.mass import calculate_mass
from Bio import SeqIO

def total_cleavage(sequence, enzyme, missed_cl):
    """ Predicts the proteolytic fragments of a protein cleaved by an enzyme.
    
    Parameters
    ----------
    sequence : str
        The sequence of the polypeptide to be cleaved.
    enzyme : str
        The enzyme that performs the cleavage. The cleavage is performed 
        following the ExPasy rules. Default is 'trypsin'.
        
    Returns
    -------
    out : dict
        A dictionary with the sequences of the cleaved peptides as keys
        and the masses of these peptides as values.
    """
    all_peptides = list(cleave(sequence, expasy_rules[enzyme], missed_cl))
    out = {seq:calculate_mass(seq,average=False) for seq in all_peptides}
    return out

def total_cleavage_from_fasta(protein_file, enzyme, missed_cl):
    """ Predicts the proteolytic fragments of a protein or a set of proteins
    contained in a fasta file cleaved by an enzyme.
    
    Parameters
    ----------
    protein_file : str
        The fasta file containing the sequence(s) of the polypeptide(s) to be cleaved.
    enzyme : str
        The enzyme that performs the cleavage. The cleavage is performed 
        following the ExPasy rules. Default is 'trypsin'.
        
    Returns
    -------
    out : dict
        A dictionary with the sequences of the cleaved peptides as keys
        and the masses of these peptides as values.
    """
    cleaved_peptides = {}
    for record in SeqIO.parse(protein_file, "fasta"):
        sequence = str(record.seq)
        cleaved_peptides.update(total_cleavage(sequence, enzyme, missed_cl))
    return sequence, cleaved_peptides

def cleavage_multiple_charge_states(fasta_file,enzyme='trypsin',charge_states=(1,5)):
    """ Predicts the proteolytic fragments of polypeptide(s) cleaved by an enzyme
    and assigns m/z values for a range of charges.The polypeptide(s) are stored 
    in a fasta file. 
    
    Parameters
    ----------
    fasta_file : str
        The fasta file containing the sequence(s) of the polypeptide(s) to be cleaved.
    enzyme : str
        The enzyme that performs the cleavage. The cleavage is performed 
        following the ExPasy rules. Default is 'trypsin'.
    max_charge_state : int
        The m/z values are assigned to each cleaved peptide for a range of charge 
        states going from z=1 to z=max_charge_state. Default is 5.
        
    Returns
    -------
    cleaved_peptide : dict
        keys: sequences of the cleaved peptides.
        values: list of masses at increasing m/z ratios
    """
    cleaved_peptides = total_cleavage_from_fasta(fasta_file,enzyme)
    for peptide in list(cleaved_peptides.keys()):
        masses = []
        for charge in range(charge_states[0],charge_states[1]+1):
            masses.append(calculate_mass(peptide,charge=charge))
        cleaved_peptides[peptide] = masses
    return cleaved_peptides


#cp = total_cleavage_from_fasta('lyso.fasta','trypsin')
#cp_df = pd.DataFrame(np.array((list(cp.keys()),list(cp.values()))).T)
#cp_df.to_excel('lyso_tryp.xlsx')
    
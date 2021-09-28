from calculate_cleavage import *
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import pandas as pd

def create_peptide_map(fasta, mc, r, min_pep):
    """
    Simulates a peptide map of a HDX-MS experiment.

    Parameters
    ----------
    fasta : str
        Fasta file containing the sequence of the protein
    mc : int
        Maximum number of missed cleavages
    r : float
        Random parameter driving the loss of peptides due to random effects.
        0 < r < 1
        r=0: every peptide is found; r=1: no peptide is found
    min_pep : int
        Minimum length of a peptide

    Returns
    -------
    None.    
    """

    out = fasta_file.split('.')[0]

    s, p = total_cleavage_from_fasta(fasta_file, 'pepsin ph2.0', missed_cl=mc)

    assignments = []; k=0
    residue_coverage = np.zeros(len(s), dtype=int)
    for pep in list(p.keys()):
        ass = []    
        if len(pep) > min_pep and random()>r:
            start = s.find(pep)+1
            end = start+len(pep)-1
            assignments.append([k,start,end,pep])
            residue_coverage[start-1:end-1] += 1
            k += 1
        
    df = pd.DataFrame(assignments)
    df.columns = ['index','start','end','seq']
    df = df.sort_values(by=['start','end'])
    df = df.reset_index(drop=True)
    df['index'] = [i+1 for i in range(df.shape[0])]

    num_of_pep = len(df)
    seq_cov = round(np.count_nonzero(residue_coverage)/len(s)*100,1)
    av_pep_len = round(np.average([len(df['seq'][i]) for i in range(df.shape[0])]),1)
    red = round(np.average(residue_coverage),1)
    text = 'Number of peptides: '+str(num_of_pep)+'\n'+\
        'Coverage: '+str(seq_cov)+'\n'+\
        'Average peptide length: '+str(av_pep_len)+'\n'+\
        'Redundancy: '+str(red)

    plt.figure()
    for i in range(df.shape[0]):
        plt.hlines(i+1,df['start'][i],df['end'][i],color='black')  
        plt.text(1,num_of_pep,text,va='top',
                 bbox={'facecolor': 'grey', 'alpha': 0.2, 'pad': 10})
    plt.xlabel('Residue Index')
    plt.ylabel('Peptide Index')
    plt.show()
    
    plt.savefig(out+'.peptidemap.png')
    df.to_csv(out+'.list',sep=' ',header=None,index=False)
    df[['index','start','end']].to_csv(out+'.ass',sep=' ',header=None,index=False)
    

if __name__ == '__main__':
   
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--fasta")
    parser.add_argument("--mc")
    parser.add_argument("--r")
    parser.add_argument("--mpl")

    opts = parser.parse_args()

    if opts.fasta:
        fasta_file = opts.fasta
    
    if opts.mc:
        mc = int(opts.mc)
    else:
        mc = 0
    
    if opts.r:
        r = float(opts.r)
    else:
        r = 0
    
    if opts.mpl:
        min_pep = int(opts.mpl)
    else:
        min_pep = 1
    
    create_peptide_map(fasta_file, mc, r, min_pep)
    
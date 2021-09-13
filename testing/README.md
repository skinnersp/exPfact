Here we report inputs to reproduce the results for the test case reported in Skinner et al, Biophys J 116, 1–10, (2019) https://doi.org/ 10.1016/j.bpj.2019.02.024


Reference protection factors are in ref.pfact 
Assignments for the sequence in test.seq are in test.list
Generate reference deuterium uptake:
% python ../python/pfact2dpred.py --ass test.list --pfact test.pfact --pH 7 --temp 300 --times test.times --seq test.seq --out test
mv test.Dpred test.Dexp
test.Dexp: experimental/reference deuterium uptake (0 <= D <=1)
format: time(hours) D(peptide1) D(peptide2) .... D(peptideN)

Predict protection factors from the reference deuterium uptake

% python ../python/exPfact.py --temp 300 --pH 7.0 --dexp test.Dexp --ass test.list --weights test.weights --tol 1e-11 --rand 10000 --seq test.seq --out out

here weights test.weights are irrelevant (all=1) and can be omitted; use the inverse of the standard deviation on measured deuterium uptake in case these are available

the file out.Dpred should be identical to test.Dexp

the file out.pfact is the set of protection factors compatible with out.Dpred (not identical to test.pfact because of the underdetermination of the problem)

% test.sh
produces 100 predictions of protection factors sets all reproducing exactly the reference deuterium uptake
variability in the prediction for specific residues provides a measure of the underdetermination of the problem. Use different assignment sets to see how the overlap of assigned peptides affect the quality of the prediction.

input: 

the output files are 
out.pfact  
format: residue_number lnP(predicted)

out.Dpred 
format: as test.Dexp with values consistent with predicted lnP

out.diff
format: peptide_id RMSD(Dpred-Dexp)
In the directory utils are some scripts that perform the calculations in the paper and the analysis.
Scripts are in bash, python, R and a bit of a mess. Simon and/or Gael will make them nice enough to put in the repository in due course. If anything is missing or confusing, please let us know: e.paci@leeds.ac.uk or s.p.skinner@leeds.ac.uk. We'd be happy to help.

% ../sh/analyse.sh test test.seq 50
analyses the results (i.e., the files out[1-100].pfact, calculating averages/median/etc of the 50% with lowest score)

Produces mean, median, interquartile range (minmax) etc. Extension cpfact includes only residues for which a pfactor is defined. I also creates a directory singlep with distribution of protection factors for individual residues.

% ../sh/do_multi.sh test
performs the multivariate analysis of the results; in the directory singlep_contig there will be files n-m.mclust that have mean and standard deviation for each of the contiguous fragments (i.e., not interrupted by gaps) going from residue n to m.

do_multi.sh uses R; you need to install a specific library by running first R interactively and 
> install.packages("mclust")

% ../sh/iso.sh
file.pfact: calculates the isotopic envelope for a given fragment using the protection factors in file.pfact
in the example peptide 2 is used. the isotopic envelope of the undeuterated sequence should be obtained from, e.g., 
http://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msisotope
need also to compile iso.f
% cd fortran ; gfortran iso.f -o iso

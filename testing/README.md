# Tutorial

Here we report the procedure to reproduce the results for the test case reported in  [Skinner et al., 'Estimating constraints for protection factors from HDX-MS data' (2019)](https://doi.org/10.1016/j.bpj.2019.02.024).

## Generation of synthetic data

Reference protection factors are in `ref.pfact`. Format: residue ln(pfact)

Assignments for the sequence in `test.seq` are in `test.list`. Format: index start end sequence

To generate reference deuterium uptake (without error):

``` python ../python/pfact2dpred.py --ass test.list --pfact test.pfact --pH 7 --temp 300 --times test.times --seq test.seq --out test ```

The script returns a file `test.Dpred` containing the predicted deuterium uptake. Format: time(hours) D(peptide1) D(peptide2) ... D(peptideN)

Please rename the file `test.Dpred` as `test.Dexp`:
* `mv test.Dpred test.Dexp` if you use Linux
* `REN test.Dpred test.Dexp` if you use Windows

## Predict protection factors from synthetic data

To find a solution, i.e. a pattern of protection factors, which accurately fits the data, ExPfact performs a random search followed by a least-squares minimization.
The random search consists in i) randomly initializing `--rand` sets of protection factors with the constraint `0<ln(Pfact)<=20`, ii) calculate the cost function (sum of squared residuals) for each set, and iii) select the set with lowest cost function.
The selected pattern of protection factors is used as initial guess for the minimization. 

To calculate one solution:

` python ../python/exPfact.py --temp 300 --pH 7.0 --dexp test.Dexp --ass test.list --weights test.weights --harm 0 --rand 10000 --seq test.seq --out out `

Parameters:
* `--temp 300`: temperature of the experiment in K 
* `--pH 7`: pH of the experiment; please insert the corrected `pH = pHread + 0.4` 
* `--dexp test.Dexp`: file containing the experimental deuterium uptake for all peptide at all time points
* `--ass test.list`: peptide assignments
* `--weights test.weights`: weights to be assigned to each measure. For synthetic data with no error, these weights are irrelevant (all set to 1). If available, use the invere of the standard deviation on measured deuterium uptake. 
* `--harm 0`: if different from 0, it introduces a penalty term to the cost function (see below)
* `--rand 10000`: number of random sets of protection factors in the random search.
* `--seq test.seq`: file containing the sequence of the protein.
* `--out out`: name used for output files

The script generates three output files:
* `out.pfact`: the predicted pattern of protection factors. 
* `out.Dpred`: the predicted deuterium uptake.
* `out.diff`: sum of squared residuals between predicted and experimental uptake for each peptide.

In the case of synthetic data, the predicted deuterium uptake contained in the file `out.Dpred` should be almost identical to the experimental values `test.Dexp`.
However, the predicted protection factors `out.pfact` could be different from the reference `test.pfact` because of underdetermination. 
The sum of squared residuals should approach zero for all peptides in the file `out.diff`.

## Producing multiple solutions

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

## Descriptive statistics

% ../sh/analyse.sh test test.seq 50
analyses the results (i.e., the files out[1-100].pfact, calculating averages/median/etc of the 50% with lowest score)

Produces mean, median, interquartile range (minmax) etc. Extension cpfact includes only residues for which a pfactor is defined. I also creates a directory singlep with distribution of protection factors for individual residues.

## Clustering algorithm

% ../sh/do_multi.sh test
performs the multivariate analysis of the results; in the directory singlep_contig there will be files n-m.mclust that have mean and standard deviation for each of the contiguous fragments (i.e., not interrupted by gaps) going from residue n to m.

do_multi.sh uses R; you need to install a specific library by running first R interactively and 
> install.packages("mclust")

## Predicting the shape of the isotopic envelope

% ../sh/iso.sh
file.pfact: calculates the isotopic envelope for a given fragment using the protection factors in file.pfact
in the example peptide 2 is used. the isotopic envelope of the undeuterated sequence should be obtained from, e.g., 
http://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msisotope
need also to compile iso.f
% cd fortran ; gfortran iso.f -o iso

## Adding a penalization term to the cost function

In order to avoid overfitting and encourage the finding of solutions with no abrupt changes between protection factors of adjacent residues (unless guided by experimental data), 
a penalty term can be added to the cost function using the parameter `--harm`.
To choose the best value for the penalty costnant (i.e., the value of `--harm`), leave-one-out cross validation can be performed:

`python ../python/cross_validation.py --dexp test.Dexp --ass test.ass --temp 300 --pH 7.0 --seq test.seq`

The script returs a file `CV.res` containing the train and test cross-validation error (CVE) calculated at different values of `--harm`. Format: harm CVE_train CVE_test. Select the `--harm` value that minimizes the total cross validation error (CVE_train+CVE_test).

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

To generate multiple solutions, just add the parameter `--rep` to the previous command.
To produce 100 predictions of protection factors, all with similar agreement with protection factors:

` python ../python/exPfact.py --temp 300 --pH 7.0 --dexp test.Dexp --ass test.list --weights test.weights --harm 0 --rand 10000 --seq test.seq --out out --rep 100`

Variability in the prediction for specific residues provides a measure of the underdetermination of the problem. 
Use different assignment sets to see how the overlap of assigned peptides affect the quality of the prediction.

Similarly to the case of a single solution, three outputs are generated, namely `out[i].pfact`, `out[i].Dpred` and `out[i].diff`, where `[i]` identifies the i-th solution.

## Descriptive statistics

To analyse the outcomes of multiple minimization calculating average, median, minimum/maximum values, digit: 

` python ../python/descriptive.py --res out --top 50 `

Parameters:
* `--res`: prefix of resutls to be analysed (equivalent to the paramter `--out` in the command to generate multiple solutions). E.g.: if outputs out1.pfact, ..., out100.pfact have been previously generated, use `--res out`
* `--top`: percentage of best solutions to be generated (default 50%)

Outputs:
* `average.pfact`: contains the average protection factor per residue with associated standard deviation. Format: residue average(lnP) std
* `median.pfact`: contains the median protection factor per residue. Format: residue median(lnP)
* `minmax.pfact`: contains the minimum and maximum value assumed by the protection factor. Format: residue min(lnP) max(lnP)
* `all.sp`: contains all the protection factors estimated by the best `--top` runs. Format:

|       | Residue 1 | Residue 2 | ... | Residue N |
| ----- | --------- | --------- | --- | --------- |
| Run 1 | lnP(1,1)  | lnP(1,2)  | ... | lnP(1,N)  |
| Run 2 | lnP(2,1)  | lnP(2,2)  | ... | lnP(2,N)  |
|  ...  |    ...    |    ...    | ... |    ...    |

Each column in the file `all.sp` can be used to build the histogram of protection factors predicted for a specific residue. 

## Clustering algorithm

The clustering algorithm is applied to the results in the file `all.sp` previously generated. 

The clustering algorithm is based on Gaussian mixture models. 
Regions covered by contiguous overlapping peptides are considered one at a time. 
The histograms of the protection factors of every residue in this region are combined into a multi-dimensional probability distribution. 
The probability distribution is fitted with a mixture of gaussians with number of components varying from 1 to 99. 
The final number of components is decided based on the Bayesian Information Criterion (BIC). 

To run the clustering algorithm: 

``` python ../python/clustering.py --ass test.ass ```

**Note**: if not done previously, run R interactively and install the mclust package: `install.packages("mclust")`

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

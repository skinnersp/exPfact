# ExPfact

ExPfact is a statistically robust approach to estimate protection factors at resolution of the single amide from sparse and underdetermined HDX-MS data. Depending on the quality of the dataset analysed, ExPfact returns one or a discrete number of possible patterns of protection factors in agreement with experimental data. 

If you use ExPfact in your work, please cite [Skinner et al., 'Estimating constraints for protection factors from HDX-MS data' (2019)](https://doi.org/10.1016/j.bpj.2019.02.024)

ExPfact is licensed to academic users under the GNU GPL Version 2 licence. Commercial users should contact Dr Emmanuele Paci: e.paci@unibo.it.

## Install

Install [Anaconda3](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).

Create a new environment using python version 3.6. Open the anaconda prompt and type:

``` conda create --name python36 python=3.6 ``` 

```  conda activate python36 ```

Install the dependencies (it is **crucial** to install the correct versions as shown below): 

``` conda install numpy=1.16.1 ```

``` conda install -c anaconda scipy=1.0.1 ```

``` conda install -c anaconda cython ```

``` conda install pandas ```

``` pip install pyopenms ```

``` conda install -c conda-forge scikit-learn ```

Before running ExPfact, you also need to compile cython code:

``` cd python ```

``` python setup_calc_dpred.py build_ext --inplace ```

The clustering algorithm is implemented in the R package [mclust](). To install R:

``` conda install -c conda-forge r-base ```

Then, to install the package, run interactively R from the command line (just digit `R`), then:

``` install.packages("mclust") ```

## Getting started

It is highly recommended to accurately follow the tutorial on synthetic data before running ExPfact on experimental data. The tutorial is contained in the **testing** repository: please follow instructions shown in the readme file. 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

dexp  = pd.read_csv("moprp.dexp",header=None,delim_whitespace=True)
weig  = pd.read_csv("moprp.weights",header=None,delim_whitespace=True)

diff  = pd.read_csv("test.diff",header=None,delim_whitespace=True)
dpred = pd.read_csv("test.Dpred",header=None,delim_whitespace=True)
pfact = pd.read_csv("test.pfact",header=None,delim_whitespace=True)

def read_dexp(dexp_file):
    """
    Reads dexp values and time points from dexp_file and returns them as numpy arrays
    :param dexp_file:
    :return:
    """
    frag_data = [line.strip().split()[1:] for line in open(dexp_file, 'r').readlines()]
    time_data = [line.strip().split()[:1] for line in open(dexp_file, 'r').readlines()]
    array = []
    for row in frag_data:
        a = [float(x) for x in row]
        array.append(a)
    time_points = np.array([float(x) for x in time_data for x in x])
    return np.array(array).T, time_points

for i in range(1,dexp.shape[1]):
    plt.figure()
    plt.errorbar(dexp[0],dexp[i],yerr=[1/weig[i][j] for j in range(len(weig[i]))],
        fmt='o',color='black',capsize=3)
    plt.plot(dpred[0],dpred[i],'--',color='black')
    plt.xscale("log")
    plt.ylim(0,1)
    plt.xlabel("Time [hr]",fontsize=12)
    plt.ylabel("Deut Uptake",fontsize=12)
    plt.title("Pep "+str(i),fontsize=12)
    plt.savefig("images/pep"+str(i)+".png")
    #plt.show()

plt.figure()
for i in range(len(diff)):
    plt.vlines(diff[0][i],ymin=0,ymax=diff[1][i],color='black')
    plt.xlabel("Peptide Index",fontsize=12)
    plt.ylabel("SSR",fontsize=12)
    plt.savefig("images/testdiff.png")
#plt.show()

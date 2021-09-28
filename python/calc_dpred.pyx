"""
Copyright (C) 2019-2020 Simon P. Skinner

This program is free software: you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from math import exp
import numpy as np

def get_rate_res(double [:] kint, double[:] P):

    cdef int N = kint.shape[0]
    cdef kint_out = np.zeros(N)
    cdef int i

    for i in range(N):
        if kint[i] == -1:
            kint_out[i] += -1
        else:
            kint_out[i] += kint[i] / exp(P[i])

    return kint_out



def res2frag(double [:] rate_res, long [:, :] assignments, double [:] time_points):

    cdef int nfrag = assignments.shape[0]
    cdef int ntime = time_points.shape[0]
    cdef double [:, :] rate_frag = np.zeros((nfrag,ntime))
    cdef double [:] namide = np.zeros(nfrag)
    cdef int i
    cdef int j

    for i in range(nfrag):
        for j in range(assignments[i][1],assignments[i][2]):
            if  (rate_res[j]>=0):
                namide[i]=namide[i]+1

    for i in range(nfrag):
        for k in range(ntime):
            for j in range(assignments[i][1],assignments[i][2]):
                if (rate_res[j]>=0):
                    rate_frag[i][k]=rate_frag[i][k]+exp(-rate_res[j]*time_points[k])
            rate_frag[i][k]=(namide[i]-rate_frag[i][k])/namide[i]

    return rate_frag


def calculate_dpred(double [:] P,
                    double [:] time_points,
                    double [:] kint,
                    long [:, :] assignments):

    rate_res = get_rate_res(kint, P)

    frag_rates = res2frag(rate_res, assignments, time_points)

    return frag_rates
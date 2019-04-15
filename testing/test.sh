#!/bin/bash

for i in `seq 1 1 100` ; do

python3 ~/local/exPfact-master_seq/python/exPfact.py --temp 300 --pH 7.0  --dexp test.Dexp --ass test.list --weights test.weights --tol 1e-11 --rand 10000 --seq test.seq --out out$i

done



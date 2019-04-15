#!/bin/bash 
# get isotopic envelopes for a given peptide e.g. with
# http://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msisotope
# set m=fragment
# set t=time point in the file t.times

t=test
m=2

ass=test.ass
kint=test.kint
awk '$2=="" {$2="0"} {print $1,$2}' $kint > tmp.kint
awk '$1!="#" && $3>0 {print $1,$3/100,$2}' $m.prospector.ucsf.edu > $m.iso
echo $m > iso.inp
echo $ass >> iso.inp
echo tmp.kint >> iso.inp
echo ${1} >> iso.inp
echo $t.times >> iso.inp
../fortran/iso < iso.inp

n=`wc out.Dpred | awk '{print $2/$1-1}'`
for j in `seq 1 1 3` ; do
    k=`awk 'NR=='$j' {for (i=2;i<=NF-1;i++){printf "%f,",$i};print $NF}' out.Dpred`
    tp=`awk 'NR=='$j' {printf "%f",$1}' out.Dpred`
    
    echo $j $n $tp 
    
    sed -e 's/KKK/'$k'/g' ../R/binomial.tmpl > r.inp
    
    R --no-save < r.inp >& r.out
    
    awk '$1=="Pr" {print $2,$3}' r.out > $j.dn
    
    cat $m.iso $j.dn >  tmp
    
    awk 'NF==2 {a[$1]=$2} NF==3 {b[$3]=$2} END {for (i in a) {for (j in b) {print i*1.00627+j,a[i]*b[j]}}}' tmp | sort -n > tmp1
    awk '{a[$1]=a[$1]+$2} END {for (i in a){print i,a[i]}}' tmp1 | sort -n > ${1}.$m.$j.ist

done
mv out.Dpred ${1}.$m.Dpred
#rm *.dn tmp1 tmp tmp.kint r.inp r.out 

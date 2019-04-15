#!/bin/bash 

awk '{printf "%d\n",$1}' costave.pfact > tmp1$$
nmax=`wc costave.pfact| awk '{print $1}' `
let "nmax = $nmax - 1"
echo $nmax
mkdir -p singlep_contig
awk '{print $2+1,$3}' ${1}.ass > tmp9$$
awk 'BEGIN {max=-1} $1>max {printf "%d \n %d ",max,$1} $2>=max {max=$2} END {print max}' tmp9$$ | awk '$1!=-1 {print $0}' > tmp0$$  

n=`awk 'END {print NR}' tmp0$$`
for l in `seq 1 1 $n` ; do
    
    ii=`awk 'NR=='$l' {print $1}' tmp0$$`
    jj=`awk 'NR=='$l' {print $2}' tmp0$$`
    i=`awk '$1=='$ii' {print NR}' tmp1$$`
    j=`awk '$1=='$jj' {print NR}' tmp1$$`

    echo $l $i $j $ii $jj

    sed -e 's/XXX/"'$i'"/g' -e 's/YYY/"'$j'"/g' ../R/multi.r > tmp.r
    R --no-save < tmp.r >& r.out 
    awk 'NR>1 {print $0}' tmp.mod > tmp2$$
    sed -e 's/\"//g' -e 's/V//g' tmp2$$ > tmpmod$$
    awk 'NR>1 {print $0}' tmp.pro > tmp3$$ 
    sed -e 's/\"//g' tmp3$$ | sort -n -k 2 -r > singlep_contig/$i"-"$j.pro
    awk 'NR>='$i' && NR<='$j' {print $0}' tmp1$$ > tmp1a$$
    paste tmp1a$$ tmpmod$$ | awk '{$2="";print $0}' > singlep_contig/$i"-"$j.mod
    awk 'NR>1 {print $0}' tmp.var > tmp2$$
    sed -e 's/\"//g' -e 's/V//g' tmp2$$ > tmpvar$$
    awk 'NR>='$i' && NR<='$j' {print $0}' tmp1$$ > tmp1a$$
    paste tmp1a$$ tmpvar$$ | awk '{$2="";print $0}' > singlep_contig/$i"-"$j.var
    mv tmp.bic singlep_contig/$i"-"$j.bic
    mv tmp.z singlep_contig/$i"-"$j.z

    cd singlep_contig
    n=`awk 'NR==2 {print NF-1}' $i"-"$j.mod`

    echo $n "components"
    rm -rf tmp$$ tmp1$$

    n=$(( $n + 1 ))
    for k in `seq 2 1 $n` ; do

    m=$(( $k -2))
    n=$(( $k -1))

#    awk 'BEGIN {print "@ s'$m' legend \"component '$n'\""} {print $1,$'$k'}' $i"-"$j.mod >> tmp$$
#    awk 'BEGIN {print " "} {print $'$k'}' $i"-"$j.var >> tmp1$$
    awk '{print $1,$'$k'}' $i"-"$j.mod >> tmp$$
    awk '{print $'$k'}' $i"-"$j.var >> tmp1$$

    echo "&" >> tmp$$
    echo "" >> tmp1$$

    done
    paste tmp$$ tmp1$$  >  $i"-"$j.mclust

    cd ..
    
done
rm tmp*$$
exit

#!/bin/bash 

if [ "${3}" = "" ]
then
   echo "need an argument .assignment/.Dexp & .seq & top xx% selected" ; exit
fi

n=`wc ${1}.ass | awk '{print $1}'`
d=`wc ${1}.Dexp | awk '{print $1}'`
HOME=`pwd`
thresh=`echo  100 / ${3} | bc -l `
echo $thresh

echo $HOME

echo $dir
for i in out*.diff ; do
    awk '{a+=$2*$2} END {print FILENAME,'$d'*a/'$n'}' $i >> tmp$$    
done
#    awk '$2<="'$thresh'" {print $0}' tmp$$ | sort -n -k 2 > tmp1$$
# change so that thresh becomes irrelevant: pick the best 10%
# check why you multiply diff by d 
     sort -g -k 2 tmp$$ > tmpz$$
     fs=`wc -l tmpz$$ | awk '{print $1}'`
     fs=`echo $fs / $thresh | bc -l `
     awk 'NR<='$fs' {print $0}' tmpz$$ > tmp1$$
     thresh=$fs
    echo -n "samples with diff less than " $thresh " "; wc tmp1$$ | awk '{printf "%s ", $1}' ; echo -n "out of " ; wc tmp$$ | awk '{print $1}'
    echo -n "best sample " ; head -1 tmp1$$
    echo -n "average cost for all " ; awk '{a+=$2;n+=1} END {print a/n}' tmp$$
    awk '{printf "%s ", $1}' tmp1$$ > tmp1a$$
    cat `cat tmp1a$$` > tmp3$$
    awk '{a[$1]+=$2;n[$1]++} END {for (i in a){print i,a[i]/n[i]}}' tmp3$$ | sort -g > diff.costave
    sed -e 's/diff/pfact/g' tmp1a$$ > tmp2$$
    cat `cat tmp2$$` > tmp3$$
    awk '{a[$1]+=$2;b[$1]+=$2*$2;n[$1]++} END {for (i in a){print i,a[i]/n[i],sqrt(b[i]/n[i]-(a[i]/n[i])^2)}}' tmp3$$ | sort -g > costave.pfact

    mkdir -p singlep 
    touch tmp4$$
    rm -f tmp1x$$ tmp2x$$
    for i in `awk '{printf "%d ",$1}' costave.pfact` ; do
	awk '$1=='$i' {printf "%9.6f \n",$2}' tmp3$$ > singlep/$i.sp
	paste tmp4$$ singlep/$i.sp > tmp5$$
	mv tmp5$$ tmp4$$
	awk 'BEGIN {for (i=1;i<=100;i++) a[i]=0} {a[int($1*5)]++;n++} END {for (i=1;i<=100;i++) {print (i+0.5)/5,a[i]*5/n}}' singlep/$i.sp | sort -g > singlep/$i.pr
        sort -g singlep/$i.sp > tmpx$$
        nmax=`wc singlep/$i.sp | awk '{print $1}' `
        awk 'BEGIN {max=-1;min=21} NR>=int("'$nmax'"*0.25) && NR<=int("'$nmax'"*0.75) {{if ($1>=max){max=$1}}; {if ($1<=min){min=$1}}} NR<=int("'$nmax'"*0.5) {med=$1} END {print "'$i'",med}' tmpx$$ >> tmp1x$$
        awk 'BEGIN {max=-1;min=21} NR>=int("'$nmax'"*0.25) && NR<=int("'$nmax'"*0.75) {{if ($1>=max){max=$1}}; {if ($1<=min){min=$1}}} END {print "'$i'",(max+min)/2,(max-min)/2}' tmpx$$ >> tmp2x$$
    done
    mv tmp4$$ singlep/all.sp
    echo  $dir " end "
    echo " "
    mv tmp$$ diff.list
    mv tmp2$$ pfact.list

    awk 'NF==2 {print $0}' tmp1x$$ | sort -g > median.pfact
    awk 'NF==3 && $3!=-11 {print $0}' tmp2x$$ | sort -g > minmax.pfact
    rm -f tmp*$$
    cd $HOME


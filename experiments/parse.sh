#!/bin/bash

datasets='
hamming6-2
johnson8-4-4
san200-0-9-2 
keller4
socfb-Texas84
sc-pkustk11
socfb-UF
socfb-B-anon
sc-pwtk
soc-flixster
tech-WHOIS
ia-wiki-Talk
socfb-Duke14
socfb-Indiana
socfb-OR
socfb-Texas84
sc-msdoor
soc-FourSquare
sc-ldoor
soc-lastfm
soc-LiveMocha
soc-orkut
soc-BlogCatalog
soc-buzznet
soc-digg
soc-flickr
'
solutions='kplex-full 	 kplex-no-seesaw    kplex-no-rule2	 kplex-no-seesaw-rule2 kPlexS'


result=results.out
rm $result


rm $result
for k in {2..5}; do
    echo "================= k=$k =================" >> $result
    echo -en "datasets\t" >>$result
    for sol in $solutions; do
        echo -en "${sol#*-}\t" >> $result
    done
    # done
    echo >>$result
    for fn in $datasets; do
    echo -en "$fn\t" >> $result
        line=()
        for sol in $solutions; do
            outfile=$sol-$k.out
            op=`grep $fn $outfile`
            if [ -z "${op}" ]; then 
                line+=('9999999999')
            else
                op=($op)
                line+=(${op[-3]})
            fi
        done
        m=`printf "%s\n" "${line[@]}" | sort -rn | tail -n1`
        echo ${line[@]}
        for i in "${line[@]}"; do 
            if [ $i == '9999999999' ]; then
                echo -en '-\t' >> $result
            elif [ $i == $m ] ; then
                echo $i | bc -l | xargs printf "*%.1f\t" >> $result
            else
                echo $i | bc -l | xargs printf "%.1f\t" >> $result
            fi 
        done
    echo >> $result
    done
done

result=kpsize.out
rm $result
for fn in $datasets; do
echo -en "$fn\t" >> $result
for k in {2..5}; do
for sol in $solutions; do
    outfile=$sol-$k.out
    # [ -e "$fn" ] || continue
    op=`grep $fn $outfile`
    if [ -z "${op}" ]; then 
        echo -en "-\t" >> $result
    else
        op=($op)
        echo -en "${op[-5]}\t" >> $result
    fi
done
done
echo >> $result
done


result=seesaw.out
rm $result
sol=kplex-full
for k in {2..5}; do
    echo >> $result

    echo "seesaw-cost">>$result
    for fn in $datasets; do

            outfile=$sol-$k.out
            # [ -e "$fn" ] || continue
            op=`grep $fn $outfile`
            if [ -z "${op}" ]; then 
                echo -en "-\t" >> $result
            else
                op=($op)
                echo -en "${op[-1]}\t" >> $result
            fi
            echo >> $result
    done
done
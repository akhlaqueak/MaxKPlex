#!/bin/bash
               
#SBATCH --job-name=kplex             ### Name of the job                
#SBATCH --ntasks=96                    ### Number of Tasks
#SBATCH --cpus-per-task=1              ### Number of Tasks per CPU    
#SBATCH --mem-per-cpu=10000                        ### Memory required, 4 gigabyte
#SBATCH --partition=medium             ### Cheaha Partition    
#SBATCH --time=24:00:00                 ### Estimated Time of Completion, 1 hour
#SBATCH --output=%x.out              ### Slurm Output file, %x is job name, %j is job id    
#SBATCH --error=%x.err               ### Slurm Error file, %x is job name, %j is job id    
#SBATCH --mail-type=END
#SBATCH --mail-user=akhlaque.ak@gmail.com

# datasets=`ls ~/graphs/network_repo/*.bin ~/graphs/dimacs10/*.bin ~/graphs/dimacs/*.bin | xargs -n 1 basename`
# datasets=`ls ~/graphs/network_repo/*.bin ~/graphs/dimacs10/*.bin ~/graphs/dimacs/*.bin`
datasets=`ls ~/graphs/network_repo/*.bin | xargs -n 1 basename`
solutions='kplex-full KpLeX kPlexS Maple kPlexT'


# solutions='kplex-full kplex-doubt kPlexS Maple'
result=results.out
rm $result

# for k in {2..5}; do
#     echo "================= k=$k =================" >> $result
#     echo -en "datasets\t" >>$result
#     for sol in $solutions; do
#         echo -en "${sol#*-}\t" >> $result
#     done
#     # done
#     echo >>$result
#     for fn in $datasets; do
#     echo -en "$fn\t" >> $result
#         line=()
#         for sol in $solutions; do
#             outfile=$sol-$k.out
#             op=`grep $fn $outfile`
#             if [ -z "${op}" ]; then 
#                 line+=('9999999999')
#             else
#                 op=($op)
#                 line+=(${op[-3]})
#             fi
#         done
#         m=`printf "%s\n" "${line[@]}" | sort -rn | tail -n1`
#         sm=`printf "%s\n" "${line[@]}" | sort -rn | tail -n2|head -n1`
#         echo ${line[@]}
#         for i in "${line[@]}"; do 
#             if [ $i == '9999999999' ]; then
#                 echo -en '-\t' >> $result
#             elif [ $i == $m ] ; then
#                 echo $i | bc -l | xargs printf "*%.1f\t" >> $result
#             elif [ $i == $sm ] ; then
#                 echo $i | bc -l | xargs printf "+%.1f\t" >> $result
#             else
#                 echo $i | bc -l | xargs printf "%.1f\t" >> $result
#             fi 
#         done
#     echo >> $result
#     done
# done

# solutions='kplex-full kplex-doubt kPlexS Maple'
result=kpsize.out
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
        for sol in $solutions; do
            outfile=$sol-$k.out
            op=`grep $fn $outfile`
            if [ -z "${op}" ]; then 
                echo -en '-\t' >> $result
            else
                op=($op)
                echo -en "${op[-5]}\t" >> $result
            fi
        done
    echo >> $result
    done
done

declare -A res
result=summary.out
ticks='1 3 10 100 300 1000 1800'
rm $result
for k in {2..5}; do
    for sol in $solutions; do
        for tick in $ticks; do
            res[$k,$sol,$tick]=0
        done
    done
done



for k in {2..5}; do
    for fn in $datasets; do
        for sol in $solutions; do
            outfile=$sol-$k.out
            op=`grep $fn $outfile`
            if ! [ -z "${op}" ]; then 
                op=($op)
                for tick in $ticks; do
                    if (( $(echo "${op[-3]} < $tick" | bc -l) )); then
                        (( res[$k,$sol,$tick]++ ))
                    fi
                done
            fi
        done
    done
done
    
echo -en "datasets\t" >> $result
for tick in $ticks; do
        echo -en "$tick\t" >> $result
done
echo >> $result

for k in {2..5}; do
    echo "================= k=$k =================" >> $result
    for sol in $solutions; do
        echo -en "$sol\t" >> $result
        for tick in $ticks; do        
            echo -en "${res[$k,$sol,$tick]}\t" >> $result
        done
        echo >> $result
    done
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

calculate_on_cutoff(){
datasets=`ls ~/graphs/$1/*.bin | xargs -n 1 basename`
result="cutoff-$1.out"
solutions='kplex-full KpLeX kPlexS Maple kPlexT'
rm $result
for k in {2..5}; do
    for sol in $solutions; do
        for tick in $ticks; do
            res[$k,$sol]=0
        done
    done
done



for k in {2..5}; do
    for fn in $datasets; do
        for sol in $solutions; do
            outfile=$sol-$k.out
            op=`grep $fn $outfile`
            if ! [ -z "${op}" ]; then 
                (( res[$k,$sol]++ ))
            fi
        done
    done
done
    
echo -en "solutions\t" >> $result
for k in {2..5}; do
        echo -en "$k\t" >> $result
done
echo >> $result

for sol in $solutions; do
    echo -en "$sol\t" >> $result   
    for k in {2..5}; do
        echo -en "${res[$k,$sol]}\t" >> $result
    done
    echo >> $result
done
cat $result
}

for dim in dimacs dimacs10 network_repo; do
    echo ==============$dim================
    echo ls ~/graphs/$dim/*.bin | wc
    # calculate_on_cutoff $dim
done
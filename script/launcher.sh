#!/bin/bash
# arg #1 : number of CPUs
# arg #2 : number of runs
# arg #3 : prefix tag i.e. _RES
# arg #4 : file name for result file 
# arg #5 : command
nb_cpus=$1
shift
nb_runs=$1
shift
tag=$1
shift
res_file=$1
shift
command=$*

function run()
{
    sleep 1
    
    msg=$1
    shift
    cpu_id=$1
    shift
    nb_runs=$1
    shift
    tag=$1
    shift
    res_file=$1
    shift
    command=$*
    
    #source /shared/softs/cuda-5.0/setenv
    
    for (( run_i=1; run_i<=${nb_runs}; run_i++))
    do
        echo ${cpu_id} ${run_i}
        ${command} | grep ${tag} | cut -d ' ' -f 2 >> "${res_file}-${cpu_id}"
        #${command} | grep "_SOR" >> "${res_file}-${cpu_id}"
    done
}


function ChildReturned()
{
    echo "A child returned"
    jobs -l
}


rm -f ${res_file}*
ts_start=$(date +%s.%N)

# Launch sub-processes
for (( cpu_i=0; cpu_i<nb_cpus; cpu_i++ ))
do
    (run "${msg}" ${cpu_i} ${nb_runs} ${tag} ${res_file} ${command}) &
    #taskset -pc ${cpu_i},$((${cpu_i}+4)) $!
    taskset -pc ${cpu_i} $!
done

# Wait for all sub-processes to finish
for i in $(jobs -p)
do
    wait $i
done

ts_stop=$(date +%s.%N)
echo "${ts_stop}-${ts_start}" | bc

# Concatenate the results
for (( cpu_i=0; cpu_i<nb_cpus; cpu_i++ ))
do
    cat ${res_file}-${cpu_i} >> ${res_file}
    rm ${res_file}-${cpu_i}
done

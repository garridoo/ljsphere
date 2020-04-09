#!/bin/bash

# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

parallel() {

    # default parameters
    local cmd_list=$1
    local n_cores=$2
    local sleep_time=2
       
    # initializations
    array[$n_cores]=pid
    
    local n=1
    local n_total=$(wc -l $cmd_list | sed 's/ .*//g')
    
    # fork first batch of processes
    for ((i = 0; i <= $n_cores - 1; i++))
    do
        if [[ $n -le $n_total ]]
        then
            $(sed -n "$n"p $cmd_list) &
            pid[$i]=$!
            echo -e "thread $i started\tpid: ${pid[$i]}"
            n=$(($n + 1))
        fi
    done
    
    # while there is still data
    # if there are free cores
    # fork new processes
    while [[ $n -le $n_total ]]
    do
        sleep $sleep_time
        clear
        echo -e "n: $n\n"
        for ((i = 0; i <= $n_cores - 1; i++))
        do
            is_running=$(ps -A | grep "^ *${pid[$i]} " | wc -l)
            if [[ $is_running -eq 0 ]]
            then
                if [[ $n -le $n_total ]]
                then
                    $(sed -n "$n"p $cmd_list) &
                    pid[$i]=$!
                    n=$(($n + 1))
                    echo -e "thread $i started\tpid: ${pid[$i]}"
                fi
            else
                echo -e "thread $i running\tpid: ${pid[$i]}"
            fi
        done
    done
    
    wait

}


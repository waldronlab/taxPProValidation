#! /usr/bin/bash

physiologies=(
    'aerophilicity'
    'gram_stain'
    'biosafety_level'
    'COGEM pathogenicity rating'
    'shape'
    'spore shape'
    'arrangement'
    'hemolysis'
)

myVar=$(which R)
myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

for i in "${physiologies[@]}"
do
    (
            echo generating data for "$i"
        if [ $myRes -eq 0 ]; then
            echo "I'm not on supermicro"
            Rscript 01_run_propagation_multistate.R "$i" "$1"
        elif [ $myRes -eq 1 ]; then
            echo "I'm on supermicro"
            /usr/bin/Rscript 01_run_propagation_multistate.R "$i" "$1"
        fi

    ) &
done

wait


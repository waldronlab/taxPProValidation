#! /usr/bin/bash

physiologies=(
    'aerophilicity'
    'gram_stain'
    'biosafety_level'
    'COGEM_pathogenicity_rating'
    'shape'
    'spore_shape'
    'arrangement'
    'hemolysis'
    'halophily_class'
    'temperature_class'
    'ph_class'
)

myVar=$(which R)
myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

for i in "${physiologies[@]}"
do
    (
            echo generating data for "$i"
        if [ $myRes -eq 0 ]; then
            echo "I'm not on supermicro"
            Rscript 02_validation_multistate.R "$i" "$1"
        elif [ $myRes -eq 1 ]; then
            echo "I'm on supermicro"
            /usr/bin/Rscript 02_validation_multistate.R "$i" "$1"
        fi

    ) &
done

wait


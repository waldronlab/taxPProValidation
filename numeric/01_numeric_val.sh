#! /usr/bin/bash

physiologies=(
    'growth_temperature'
    'coding_genes'
    'genome_size'
    'length'
    'mutation_rate_per_site_per_generation'
    'mutation_rate_per_site_per_year'
    'optimal_ph'
    'width'
    'gc_content'
)

myVar=$(which R)
myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

for i in "${physiologies[@]}"
do
    (
            echo generating data for "$i"
        if [ $myRes -eq 0 ]; then
            echo "I'm not on supermicro"
            Rscript 01_numeric_val.R "$i" "$1"
        elif [ $myRes -eq 1 ]; then
            echo "I'm on supermicro"
            /usr/bin/Rscript 01_numeric_val.R "$i" "$1"
        fi

    ) &
done

wait


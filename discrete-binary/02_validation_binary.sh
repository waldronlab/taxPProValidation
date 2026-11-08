#! /usr/bin/bash

physiologies=(
   "animal_pathogen"
   "antimicrobial_sensitivity"
   "biofilm_formation"
   "extreme_environment"
   "health_associated"
   "hydrogen_gas_producing"
   "lactate_producing"
   "motility"
   "plant_pathogenicity"
   "spore_formation"
   "host-associated"
   "sphingolipid_producing"
   "butyrate_producing"
)

myVar=$(which R)
myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

for i in "${physiologies[@]}"
do
    (
            echo generating data for "$i"
        if [ $myRes -eq 0 ]; then
            echo "I'm not on supermicro"
            Rscript 02_validation_binary.R "$i" "$1"
        elif [ $myRes -eq 1 ]; then
            echo "I'm on supermicro"
            /usr/bin/Rscript 02_validation_binary.R "$i" "$1"
        fi

    ) &
done

wait


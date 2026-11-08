#!/bin/bash
# Full cross-validation: three attribute types in parallel, ranks in sequence.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=2
for d in discrete-binary discrete-multistate-intersection numeric; do
  (
    cd $d
    for rank in all genus species strain; do
      echo "== $d $rank $(date)" >> ../run_all.log
      case $d in
        discrete-binary) bash 01_prediction_binary.sh $rank > log_01_$rank.txt 2>&1 ;;
        discrete-multistate-intersection) bash 01_prediction_multistate.sh $rank > log_01_$rank.txt 2>&1 ;;
        numeric) bash 01_numeric_val.sh $rank > log_01_$rank.txt 2>&1 ;;
      esac
    done
    echo "== $d prediction done $(date)" >> ../run_all.log
  ) &
done
wait
echo "== all predictions done $(date)" >> run_all.log

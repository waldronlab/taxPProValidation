
One-liners for performing prediction and validation of numeric attributes:


```
tmux new -d -s numeric "./01_numeric_val.sh all && ./01_numeric_val.sh genus && ./01_numeric_val.sh species && ./01_numeric_val.sh strain"
/usr/bin/Rscript --vanilla 02_summary_numeric.R # path to Rscript might vary
```

Output is in the [numeric_castor-ltp_summary.tsv](./numeric_castor-ltp_summary.tsv) file.


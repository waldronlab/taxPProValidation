## Thresholds that turn scores into yes/no calls in the 02 validation scripts.
## Sourced from discrete-binary/ and discrete-multistate-intersection/ as
## source("../validation_calls.R"). Tested by test_validation_calls.R.
##
## A predicted state counts as positive at the score the export ships: 0.5 or
## more (taxPPro::scores2Freq() "sometimes" and up; ADR 0002). The held-out
## truth counts as positive only when the curated score is above 0.8, so the
## two cut-offs differ on purpose.
PRED_CALL <- 0.5
TRUTH_CALL <- 0.8
callPredicted <- function(score) !is.na(score) & score >= PRED_CALL
callTruth <- function(score) !is.na(score) & score > TRUTH_CALL
## Binary attributes: one row per taxon survives. Among rows with equal score
## the TRUE row wins, because the export ships TRUE at exactly 0.5.
pickBinaryCall <- function(df) {
    df |>
        dplyr::arrange(NCBI_ID, dplyr::desc(grepl("TRUE$", Attribute))) |>
        dplyr::group_by(NCBI_ID) |>
        dplyr::slice_max(order_by = Score, n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
}

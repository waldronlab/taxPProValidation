## Scoring functions for create_attribute_summary.R (bugphyzzExports ADR 0006).
## Tested by test_attribute_scores.R.

nsti_breaks <- c(-Inf, 0.02, 0.05, 0.1, 0.2, Inf)

## Weights matching held-out NSTI to target NSTI: share of target tips in the
## held-out tip's bin / share of held-out tips in that bin. Pass the NSTI of
## all folds together so the shares are pooled.
nsti_weights <- function(heldout_nsti, target_nsti) {
    band_h <- cut(heldout_nsti, nsti_breaks)
    share_t <- prop.table(table(cut(target_nsti, nsti_breaks)))
    share_h <- prop.table(table(band_h))
    as.numeric(share_t[as.character(band_h)] / share_h[as.character(band_h)])
}

## Cohen's kappa from weighted counts. `states` is the full set of states;
## with `positions` (same order) agreement is linear in the distance between
## positions, 1 - |p_i - p_j| / (max - min), otherwise exact agreement only.
## NA, with a warning, when agreement expected by chance is 1.
weighted_kappa <- function(pred, truth, w, states, positions = NULL) {
    stopifnot(all(pred %in% states), all(truth %in% states))
    cm <- tapply(w, list(factor(pred, states), factor(truth, states)), sum)
    cm[is.na(cm)] <- 0
    cm <- cm / sum(cm)
    agree <- if (is.null(positions)) {
        diag(length(states))
    } else {
        1 - abs(outer(positions, positions, "-")) / diff(range(positions))
    }
    expected <- sum(agree * outer(rowSums(cm), colSums(cm)))
    if (isTRUE(all.equal(expected, 1))) {
        warning("kappa undefined: chance agreement is 1", call. = FALSE)
        return(NA_real_)
    }
    (sum(agree * cm) - expected) / (1 - expected)
}

## Squared correlation, as 01_numeric_val.R reports R2, with weights.
weighted_r2 <- function(actual, predicted, w) {
    summary(lm(actual ~ predicted, weights = w))$r.squared
}

## Lower bound of a two-sided 95% t-interval over per-fold scores.
lower_bound <- function(x) {
    x <- x[!is.na(x)]
    mean(x) - qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
}

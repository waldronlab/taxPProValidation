## Attribute-level validation scores for bugphyzzExports ADR 0006.
##
## Reads the fold files the prediction scripts write (discrete: *_test_FoldN.csv
## and *_predicted_FoldN.csv; numeric: *_castor-ltp_predictions.csv) and
## writes, from the repository root:
##   validation_attributes.tsv  one row per rank and physiology
##   validation_states.tsv      per-state precision of discrete predictions
##
## Score: Cohen's kappa for discrete attributes, linearly weighted when the
## attribute's states have an ordinal_position in bugphyzz attributes.tsv;
## R2 (squared correlation, as in 01_numeric_val.R) for numeric attributes.
## Held-out tips are weighted so that their distribution of distance to the
## nearest annotated tip (NSTI) matches that of the LTP tips with no
## annotation: weight = share of target tips in the tip's NSTI bin / share of
## held-out tips in that bin, bins at 0.02, 0.05, 0.1 and 0.2, shares pooled
## over folds, no clipping. score_lcb is the lower bound of a two-sided 95%
## t-interval over the per-fold scores (k - 1 df).
##
## Usage (on supermicro, after run_all.sh):
##   Rscript create_attribute_summary.R [run_dir] [attributes.tsv]
## run_dir holds the fold files (default "."); attributes.tsv defaults to the
## installed bugphyzz's. Output is written to the current directory.

suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(castor)
    library(taxPPro)
})

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) >= 1) args[[1]] else "."
attrs_file <- if (length(args) >= 2) args[[2]] else
    system.file("extdata", "attributes.tsv", package = "bugphyzz")

nsti_breaks <- c(-Inf, 0.02, 0.05, 0.1, 0.2, Inf)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data

attrs_tsv <- read.delim(attrs_file)
stopifnot("ordinal_position" %in% names(attrs_tsv))

## Distance from each tip in `tips` to the nearest tip in `reference`.
nearest <- function(tips, reference) {
    nn <- find_nearest_tips(tree, target_tips = match(reference, tree$tip.label))
    nn$nearest_distance_per_tip[match(tips, tree$tip.label)]
}

## Weights matching held-out NSTI to target NSTI, pooled over folds.
nsti_weights <- function(heldout_nsti, target_nsti) {
    band_h <- cut(heldout_nsti, nsti_breaks)
    share_t <- prop.table(table(cut(target_nsti, nsti_breaks)))
    share_h <- prop.table(table(band_h))
    as.numeric(share_t[as.character(band_h)] / share_h[as.character(band_h)])
}

weighted_kappa <- function(pred, truth, w, states, ordinal) {
    k <- length(states)
    cm <- tapply(w, list(factor(pred, states), factor(truth, states)), sum)
    cm[is.na(cm)] <- 0
    cm <- cm / sum(cm)
    agree <- if (ordinal) {
        outer(seq_len(k), seq_len(k), function(i, j) 1 - abs(i - j) / (k - 1))
    } else {
        diag(k)
    }
    expected <- outer(rowSums(cm), colSums(cm))
    (sum(agree * cm) - sum(agree * expected)) / (1 - sum(agree * expected))
}

weighted_r2 <- function(actual, predicted, w) {
    summary(lm(actual ~ predicted, weights = w))$r.squared
}

lower_bound <- function(x) {
    mean(x) - qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
}

summarise_scores <- function(folds, base) {
    tibble(
        base,
        n_folds = nrow(folds),
        score_mean = mean(folds$score),
        score_sd = sd(folds$score),
        score_lcb = lower_bound(folds$score),
        unweighted_mean = mean(folds$unweighted)
    )
}

## Discrete ####
discrete_one <- function(dir, phys, rank) {
    stem <- file.path(run_dir, dir, paste0(phys, "_", rank, "_phytools-ltp"))
    folds <- seq_len(10)
    test_files <- paste0(stem, "_test_Fold", folds, ".csv")
    pred_files <- paste0(stem, "_predicted_Fold", folds, ".csv")
    if (!all(file.exists(c(test_files, pred_files)))) return(NULL)

    test <- map2_dfr(test_files, folds, ~ mutate(read.csv(.x), fold = .y))
    ## curated state: highest Score; ties go to the first state in name order
    truth <- test |>
        arrange(NCBI_ID, desc(Score), Attribute) |>
        distinct(NCBI_ID, .keep_all = TRUE) |>
        select(NCBI_ID, fold, truth = Attribute)
    tips <- tip_data |>
        select(tip_label, NCBI_ID) |>
        inner_join(truth, by = "NCBI_ID")
    preds <- map2_dfr(pred_files, folds, ~ {
        read.csv(.x) |>
            filter(tip_label %in% tips$tip_label[tips$fold == .y]) |>
            arrange(tip_label, desc(Score), Attribute) |>
            distinct(tip_label, .keep_all = TRUE) |>
            transmute(tip_label, pred = Attribute, fold = .y)
    })
    held <- inner_join(tips, preds, by = c("tip_label", "fold"))

    held$nsti <- NA_real_
    for (f in folds) {
        i <- held$fold == f
        held$nsti[i] <- nearest(
            held$tip_label[i], unique(tips$tip_label[tips$fold != f])
        )
    }
    targets <- setdiff(tree$tip.label, tips$tip_label)
    held$w <- nsti_weights(held$nsti, nearest(targets, unique(tips$tip_label)))

    group <- unique(test$Attribute_group)
    states <- sort(unique(c(held$truth, held$pred)))
    pos <- attrs_tsv$ordinal_position[
        match(paste(group, states), paste(attrs_tsv$attribute_group, attrs_tsv$attribute))
    ]
    ordinal <- length(states) > 2 && all(!is.na(pos))
    if (ordinal) states <- states[order(pos)]

    fold_scores <- held |>
        group_by(fold) |>
        summarise(
            score = weighted_kappa(pred, truth, w, states, ordinal),
            unweighted = weighted_kappa(pred, truth, rep(1, n()), states, ordinal)
        )
    attr_row <- summarise_scores(fold_scores, tibble(
        method = "phytools-ltp", rank = rank, physiology = group,
        metric = if (ordinal) "weighted kappa" else "kappa"
    )) |>
        mutate(
            n_heldout = nrow(held),
            ess_fraction = sum(held$w)^2 / sum(held$w^2) / nrow(held),
            max_weight = max(held$w),
            nsti_median_heldout = median(held$nsti),
            nsti_median_targets = median(nearest(targets, unique(tips$tip_label)))
        )
    state_rows <- held |>
        group_by(state = pred) |>
        summarise(
            n_predicted = n(),
            precision = sum(w * (pred == truth)) / sum(w),
            .groups = "drop"
        ) |>
        mutate(
            method = "phytools-ltp", rank = rank, physiology = group,
            state = sub("^.*--", "", state), .before = 1
        )
    list(attr = attr_row, states = state_rows)
}

## Numeric ####
numeric_one <- function(phys, rank) {
    f <- file.path(run_dir, "numeric", paste0(phys, "_", rank, "_castor-ltp_predictions.csv"))
    if (!file.exists(f)) return(NULL)
    held <- read.csv(f)
    annotated <- unique(held$tip_label)
    held$nsti <- NA_real_
    for (k in unique(held$fold)) {
        i <- held$fold == k
        held$nsti[i] <- nearest(held$tip_label[i], setdiff(annotated, held$tip_label[i]))
    }
    targets <- setdiff(tree$tip.label, annotated)
    target_nsti <- nearest(targets, annotated)
    held$w <- nsti_weights(held$nsti, target_nsti)
    fold_scores <- held |>
        group_by(fold) |>
        summarise(
            score = weighted_r2(actual, predicted, w),
            unweighted = weighted_r2(actual, predicted, rep(1, n()))
        )
    summarise_scores(fold_scores, tibble(
        method = "castor-ltp", rank = rank, physiology = gsub("_", " ", phys),
        metric = "r2"
    )) |>
        mutate(
            n_heldout = nrow(held),
            ess_fraction = sum(held$w)^2 / sum(held$w^2) / nrow(held),
            max_weight = max(held$w),
            nsti_median_heldout = median(held$nsti),
            nsti_median_targets = median(target_nsti)
        )
}

ranks <- c("all", "genus", "species", "strain")
discrete_dirs <- c("discrete-binary", "discrete-multistate-intersection")
discrete_phys <- unique(unlist(map(discrete_dirs, ~ sub(
    "_(all|genus|species|strain)_phytools-ltp_test_Fold1\\.csv$", "",
    list.files(file.path(run_dir, .x), pattern = "_test_Fold1\\.csv$")
))))
discrete <- list()
for (d in discrete_dirs) for (p in discrete_phys) for (r in ranks) {
    res <- discrete_one(d, p, r)
    if (!is.null(res)) discrete[[paste(d, p, r)]] <- res
}
numeric_phys <- unique(sub(
    "_(all|genus|species|strain)_castor-ltp_predictions\\.csv$", "",
    list.files(file.path(run_dir, "numeric"), pattern = "_predictions\\.csv$")
))
numeric <- compact(unlist(
    map(numeric_phys, function(p) map(ranks, ~ numeric_one(p, .x))),
    recursive = FALSE
))

attributes_out <- bind_rows(map(discrete, "attr"), numeric) |>
    arrange(physiology, rank) |>
    mutate(across(where(is.double), ~ round(.x, 3)))
states_out <- bind_rows(map(discrete, "states")) |>
    arrange(physiology, rank, state) |>
    mutate(precision = round(precision, 3))

stopifnot(!anyDuplicated(attributes_out[c("rank", "physiology")]))
write.table(attributes_out, "validation_attributes.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(states_out, "validation_states.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Wrote validation_attributes.tsv (", nrow(attributes_out),
        " rows) and validation_states.tsv (", nrow(states_out), " rows)")

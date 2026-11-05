
args <- commandArgs(trailingOnly = TRUE)
# args <- list('growth_temperature', 'all')

library(taxPPro)
library(dplyr)
library(bugphyzz)
library(phytools)
library(castor)
library(cvTools)
library(purrr)
library(tidyr)
library(ggplot2)
library(logr)

phys_name <- args[[1]]
rank_arg <- args[[2]]

if (rank_arg == 'all') {
    rank_var <- c('genus', 'species', 'strain')
} else {
    rank_var <- rank_arg
}

logfile <- paste0(phys_name, '_', rank_arg, '_castor-ltp')
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

msg <- paste0('Running validation for: ', phys_name)
log_print(msg)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data

## bugphyzz data ####
phys <- physiologies(gsub('_', ' ', phys_name))[[1]]

modifyNumeric <- function(x) {
    x |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown') |>
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max)) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_min))) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_max))) |>
        dplyr::mutate(Rank = taxizedb::taxid2rank(NCBI_ID, db = 'ncbi')) |> 
        dplyr::mutate(Taxon_name = taxizedb::taxid2name(NCBI_ID, db = 'ncbi')) |> 
        dplyr::mutate(NCBI_ID = addRankPrefix(NCBI_ID, Rank)) |> 
        dplyr::filter(!is.na(Rank) & !is.na(Taxon_name) & !is.na(NCBI_ID)) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
           .data$Confidence_in_curation, n = 1, with_ties = FALSE
        ) |>
        dplyr::mutate(
            Attribute_value = mean(.data$Attribute_value_min, .data$Attribute_value_max)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-Attribute_value_min, -Attribute_value_max) |>
        dplyr::distinct()
}

# dat <- modifyNumeric(phys) |> 
#     filter(Rank %in% rank_var)

dat <- phys |> 
    filterData() |> 
    getDataReady() |> 
    filter(Rank %in% rank_var)

if (!nrow(dat)) {
    log_print('not enough data for propagation')
    quit(save = 'no')
    log_close()
}
    
fdat <- dat |> 
    filter(Rank %in% rank_var) |>
    filter(NCBI_ID %in% unique(tip_data$NCBI_ID))

if (!nrow(fdat)) {
    log_print('not enough data for propagation')
    quit(save = 'no')
    log_close()
}

tip_data_annotated <- left_join(tip_data, fdat, by = c('NCBI_ID')) |>
    select(tip_label, Attribute_value) # I don't need to filter NAs here

known_tips <- tip_data_annotated |> 
    filter(!is.na(Attribute_value)) |> 
    pull(tip_label)

#ltp_per <- round(length(known_tips) / Ntip(tree) * 100)
#if (ltp_per < 1) {
#    msg <- paste0("Not enough data for ASR-LTP for: ", phys_name )
#    log_print(msg)
#    quit(save = "no")
#}

nodes_with_taxid <- grep("^n", tree$node.label, value = TRUE, invert = TRUE) ## Changed code
nsti <- getNsti(tree, known_tips, nodes_with_taxid) ## Changed code

# dist <- castor::get_all_pairwise_distances(tree = tree, only_clades = known_tips)
# dist <- dist[lower.tri(dist)]

counts <- data.frame(
    physiology = gsub('_', ' ', phys_name),
    rank = rank_arg,
    ltp_bp_phys = length(unique(fdat$NCBI_ID)),
    bp_phys = length(unique(dat$NCBI_ID)),
    ltp = length(unique(tip_data$NCBI_ID)),
    nsti_mean = round(mean(nsti$nsti), 3),
    nsti_sd = round(sd(nsti$nsti), 3)
    # dist_mean = round(mean(dist), 3),
    # dist_sd = round(sd(dist), 3)
)

total_per <- floor(length(unique(fdat$NCBI_ID)) / Ntip(tree) * 100)
if (total_per < 1) {
    message("Not enough percentage data for propagation.")
    quit(save = "no")
}

allTips <- tip_data_annotated$Attribute_value
names(allTips) <- tip_data_annotated$tip_label

values <- allTips[!is.na(allTips)]

## Create train and test sets ####
## Folds are assigned per species (strains and subspecies take their
## species' fold) and per NCBI id, not per tip: several tips can share an id
## and would otherwise put the same value in both sets (taxPProValidation#1).
set.seed(1234)
tip_keys <- tip_data |>
    select(tip_label, NCBI_ID, species_taxid) |>
    mutate(
        fold_key = ifelse(
            is.na(species_taxid) | species_taxid == '', NCBI_ID, species_taxid
        )
    )
annotated_keys <- unique(tip_keys$fold_key[tip_keys$tip_label %in% names(values)])
foldN <- cvFolds(n = length(annotated_keys), K = 10)
key_fold <- setNames(foldN$which, annotated_keys[foldN$subsets[, 1]])
tip_fold <- key_fold[tip_keys$fold_key[match(names(values), tip_keys$tip_label)]]
testSets <- vector('list', 10)
trainSets <- vector('list', 10)
for (i in 1:10) {
    foldName <- paste0('Fold', i)
    testSets[[i]] <- values[tip_fold == i]
    names(testSets)[i] <- foldName
    trainSets[[i]] <- values[tip_fold != i]
    names(trainSets)[i] <- foldName
}
stopifnot(all(map2_lgl(testSets, trainSets, ~ !any(
    tip_keys$fold_key[match(names(.x), tip_keys$tip_label)] %in%
        tip_keys$fold_key[match(names(.y), tip_keys$tip_label)]
))))

## hidden-state-prediction ####
hsp <- vector('list', 10)
for (i in seq_along(trainSets)) {
    names(hsp)[i] <- names(trainSets)[i]

    x <- testSets[[i]]
    x[!is.na(x)] <- NA
    y <- trainSets[[i]]
    z = c(x, y)
    a = allTips[!names(allTips) %in% names(z)]
    input_vector <- c(z, a)
    input_vector <- input_vector[tree$tip.label]

    res <- hsp_squared_change_parsimony(
        tree = tree, tip_states = input_vector, weighted = TRUE,
        check_input = TRUE
    )
    
    # res <- hsp_independent_contrasts(
    #     tree = tree, tip_states = input_vector, weighted = TRUE,
    #     check_input = TRUE
    # )
    
    # res <- hsp_subtree_averaging(
    #     tree = tree, tip_states = input_vector,
    #     check_input = TRUE
    # )


    statesDF <- data.frame(
        label = tree$tip.label,
        value = res$states[1:Ntip(tree)]
    )

    states <- statesDF$value
    names(states) <- statesDF$label

    hsp[[i]] <- states[names(testSets[[i]])]
}

## Calculate metrics
metrics <- map2(hsp, testSets, ~ {
    predicted_values <- .x
    actual_values <- .y

    mse <- mean((predicted_values - actual_values)^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(predicted_values - actual_values))

    ss_total <- sum((actual_values - mean(actual_values))^2)
    ss_residual <- sum((actual_values - predicted_values)^2)
    # r_squared <- 1 - (ss_residual / ss_total)
    r_squared <- summary(lm(formula = actual_values ~ predicted_values))$r.squared

    mape <- mean(abs((actual_values - predicted_values) / actual_values)) * 100

    evaluation_results <- data.frame(
        # Metric = c("MSE", "RMSE", "MAE", "R_squared", "MAPE"),
        # Metric = c("MSE", "RMSE", "r_squared"),
        Metric = c("r_squared"),
        Value = c(r_squared)
        # Value = c(mse, rmse, r_squared)
        # Value = c(mse, rmse, mae, r_squared, mape)
    )

    return(evaluation_results)
}) |>
    bind_rows() |>
    group_by(Metric) |>
    summarise(
        r_squared_mean = mean(Value, na.rm = TRUE),
        r_squared_sd = sd(Value, na.rm = TRUE)
    ) |> 
    select(-Metric) |> 
    as.data.frame()


output_tsv <- cbind(metrics, counts) |> 
    relocate(physiology) |> 
    mutate(method = 'castor-ltp')

write.table(
    x = output_tsv, file = paste0(logfile, '.tsv'), sep = '\t', quote = FALSE,
    row.names = FALSE
)

log_close()

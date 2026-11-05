
args <- commandArgs(trailingOnly = TRUE)
# args <- list('aerophilicity', 'all')
# args <- list('butyrate producing', 'all')
suppressMessages({
    library(phytools)
    library(castor)
    library(cvTools)
    library(taxPPro)
    library(bugphyzz)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(mltools)
    library(logr)
})


## Fold assignment by species (shared by the three prediction scripts).
## tip_data has species_taxid for every tip; a strain's key is its species,
## a species' or genus' key is itself.
assignSpeciesFolds <- function(dat, tip_data, K = 10) {
    keys <- tip_data |>
        mutate(NCBI_ID = as.character(NCBI_ID), species_taxid = as.character(species_taxid)) |>
        distinct(NCBI_ID, species_taxid) |>
        mutate(
            fold_key = ifelse(
                is.na(species_taxid) | species_taxid == '', NCBI_ID, species_taxid
            )
        ) |>
        distinct(NCBI_ID, .keep_all = TRUE) |>
        select(NCBI_ID, fold_key)
    dat$NCBI_ID <- as.character(dat$NCBI_ID)
    dat <- left_join(dat, keys, by = 'NCBI_ID')
    dat$fold_key[is.na(dat$fold_key)] <- dat$NCBI_ID[is.na(dat$fold_key)]
    ## one fold per key, stratified by the key's first attribute value in
    ## sorted order (deterministic for mixed keys)
    key_tbl <- dat |>
        group_by(fold_key) |>
        summarise(
            Attribute = sort(Attribute)[1], n_values = n_distinct(Attribute),
            .groups = 'drop'
        )
    if (any(key_tbl$n_values > 1)) {
        message(sum(key_tbl$n_values > 1), ' species carry more than one value ',
                'across their strains; stratified by the first in sorted order')
    }
    ## A class with fewer than K species cannot reach every fold: drop it, as
    ## remove_ids does for classes with fewer than 10 taxa. Returns NULL when
    ## fewer than two classes remain, so the caller can quit cleanly.
    small <- key_tbl |> count(Attribute) |> filter(n < K) |> pull(Attribute)
    if (length(small)) {
        message('Dropping classes with fewer than ', K, ' species: ',
                paste(small, collapse = ', '))
        key_tbl <- filter(key_tbl, !Attribute %in% small)
        dat <- filter(dat, fold_key %in% key_tbl$fold_key)
    }
    if (n_distinct(key_tbl$Attribute) < 2) return(NULL)
    folds <- split(key_tbl, key_tbl$Attribute) |>
        map(~ {
            cv <- cvFolds(n = nrow(.x), K = K)
            tibble(fold_key = .x$fold_key[cv$subsets[, 1]], fold = cv$which)
        }) |>
        bind_rows()
    out <- left_join(dat, folds, by = 'fold_key')
    stopifnot(!anyNA(out$fold))
    out
}

seed <- 42342
phys_name <- args[[1]]
rank_arg <- args[[2]]

if (rank_arg == 'all') {
    rank_var <- c('genus', 'species', 'strain')
} else {
    rank_var <- rank_arg
}

logfile <- paste0(phys_name, '_', rank_arg, '_phytools-ltp')
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

msg <- paste0('Predicting ', phys_name, '. Rank: ', rank_arg)
log_print(msg, blank_after = TRUE)

msg <- paste0('Strating at ', Sys.time())
log_print(msg, blank_after = TRUE)

## Tree data (taxPPro) ####
ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data

## bugphyzz data ####
phys <- physiologies(gsub('_', ' ', phys_name))[[1]]
dat <- getDataReady(filterData(phys)) |> 
    filter(!is.na(Evidence)) |>  
    filter(Rank %in% rank_var) 

fdat <- dat |> 
    filter(NCBI_ID %in% unique(tip_data$NCBI_ID)) |> 
    filter(!duplicated(NCBI_ID))

if (!nrow(fdat)) {
    msg <- paste0(
        'Not enough data for propagation of: ', phys_name, '; rank: ', rank_arg,
        'Quitting.'
    )
    log_print(msg)
    quit(save = 'no')
}

remove_ids <- split(fdat, fdat$Attribute) |> 
    keep(~ nrow(.x) < 10) |> 
    map(~ pull(.x, NCBI_ID)) |> 
    unlist() |> 
    unique()
fdat <- filter(fdat, !NCBI_ID %in% remove_ids)

if (nrow(fdat) < 10) {
    msg <- paste0(
        'Not enough data for propagation of: ', phys_name, '; rank: ', rank_arg,
        'Quitting. Number of attributes: ', length(keep_this)
    )
    log_print(msg)
    quit(save = 'no')
    
}

## Create folds ####
keep_this <- split(fdat, fdat$Attribute)

if (length(keep_this) < 2) {
    msg <- paste0(
        'Not enough data for propagation of: ', phys_name, '; rank: ', rank_arg,
        'Quitting. Number of attributes: ', length(keep_this)
    )
    log_print(msg)
    quit(save = 'no')
}

## Folds are assigned per species, not per taxon: strains and subspecies get
## the fold of their species so a taxon is never predicted from its own
## strains (taxPProValidation#1). Stratified by attribute value through the
## attribute of each species' first annotation.
set.seed(seed)
fdat <- assignSpeciesFolds(fdat, tip_data, K = 10)
if (is.null(fdat)) {
    msg <- paste0(
        'Not enough data for propagation of: ', phys_name, '; rank: ', rank_arg,
        '. Fewer than two classes with 10 or more species. Quitting.'
    )
    log_print(msg)
    log_close()
    quit(save = 'no')
}

test_folds <- vector('list', 10)
train_folds <- vector('list', 10)
for (i in 1:10) {
    fold_name <- paste0('Fold', i)
    test_folds[[i]] <- fdat[fdat$fold == i, ]
    names(test_folds)[i] <- fold_name
    train_folds[[i]] <- fdat[fdat$fold != i, ]
    names(train_folds)[i] <- fold_name
}

## No taxon, and no species of a tested strain, may be in both sets
stopifnot(all(map2_lgl(test_folds, train_folds, ~ !any(.x$NCBI_ID %in% .y$NCBI_ID))))
stopifnot(all(map2_lgl(test_folds, train_folds, ~ !any(.x$fold_key %in% .y$fold_key))))

## 
ltp_bp <- map_int(split(dat, dat$Attribute), ~ {
    sum(unique(.x$NCBI_ID) %in% unique(tip_data$NCBI_ID))
})

bp <- map_int(split(dat, dat$Attribute), ~ {
    length(unique(.x$NCBI_ID))
})

ltp_bp_phys <- sum(unique(tip_data$NCBI_ID) %in% unique(dat$NCBI_ID))

bp_phys <- length(unique(dat$NCBI_ID))


## Create input matrix for tree ####
known_priors <- vector('list', length(train_folds))
for (i in seq_along(train_folds)) {
    names(known_priors)[i] <- paste0('Fold', i)
    kp <- tip_data |> 
        left_join(
            train_folds[[i]], by = 'NCBI_ID', relationship = 'many-to-many'
        ) |> 
        filter(!is.na(Attribute)) |> 
        select(tip_label, Attribute, Score) |> 
        pivot_wider(
            names_from = 'Attribute', values_from = 'Score', values_fill = 0
        ) |> 
        tibble::column_to_rownames(var = 'tip_label') |> 
        as.matrix()
    
#    ltp_per <- floor(nrow(kp) / Ntip(tree) * 100)
#    if (ltp_per < 1) {
#        msg <- message("Not enoght data for ASR-LTP for: ", phys_name)
#        quit(save = "no")
#    }
    known_priors[[i]] <- kp
}

## Code for creating and exporting counts ####

nsti_input <- map(known_priors, rownames) |>
    unlist(use.names = FALSE) |>
    unique()

nodes_with_taxid <- grep("^n", tree$node.label, value = TRUE, invert = TRUE)
nsti <- getNsti(tree, nsti_input, nodes_with_taxid)$nsti

counts <- data.frame(
    Attribute = names(ltp_bp),
    ltp_bp = unname(ltp_bp),
    bp = unname(bp),
    ltp_bp_phys = ltp_bp_phys ,
    bp_phys = bp_phys,
    ltp = Ntip(tree),
    nsti_mean = round(mean(nsti), 3),
    nsti_sd = round(sd(nsti), 3)
)


## total per represents the total of the tree not each fold.
total_per <- floor(ltp_bp_phys / Ntip(tree) * 100)
if (total_per < 1) {
    message("Not enought data for propagation.")
    quit(save = "no")
}

write.table(
    x = counts, file = paste0(gsub(' ', '_', phys_name), '_', rank_arg, '_phytools-ltp_counts', '.tsv'),
    row.names = FALSE, quote = FALSE, sep = '\t'
)
####

input_mats <- vector('list', length(known_priors))
for (i in seq_along(known_priors)) {
   tips <- tree$tip.label[!tree$tip.label %in% rownames(known_priors[[i]])] 
   input_mat <- matrix(
       data = rep(1/ncol(known_priors[[i]]), length(tips) * ncol(known_priors[[i]])),
       ncol = ncol(known_priors[[i]]),
       dimnames = list(rownames = tips, colnames = colnames(known_priors[[i]]))
   )
   input_mat <- rbind(input_mat, known_priors[[i]])
   input_mat[tree$tip.label]
   names(input_mats)[i] <- paste0('Fold', i)
   input_mats[[i]] <- input_mat
}

## Run prediction/propagation ####
predictions <- vector('list', length(input_mats))
for (i in seq_along(input_mats)) {
    fold_name <- paste0('Fold', i)
    message('Running predictions for ', fold_name)
    tim <- system.time({
        fit <- fitMk(
            tree = tree, x = input_mats[[i]], model = 'ER', pi = 'equal'
        )
        names(predictions)[[i]] <- fold_name
        
        
        predictions[[i]] <- ancr(object = fit, tips = TRUE)$ace
    })
    log_print(tim, blank_after = FALSE)
}

## ASR only predictions
predictions_tips <- vector('list', length(predictions))
for (i in seq_along(predictions)) {
    x <- predictions[[i]]
    x <- x[1:Ntip(tree),]
    names(predictions_tips)[[i]] <- paste0('Fold', i)
    predictions_tips[[i]] <- x |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = 'tip_label') |> 
        pivot_longer(
        names_to = 'Attribute', values_to = 'Score', cols = 2:last_col()
        ) |> 
        left_join(tip_data, by = 'tip_label')
}


msg <- paste0('Writing files')
log_print(msg, blank_after = TRUE)

## Write test sets
for (i in seq_along(test_folds)) {
    fname <- paste0(
        gsub(' ', '_', phys_name), '_', rank_arg, '_phytools-ltp_test_Fold', i, '.csv'
    )
    write.csv(
        x = test_folds[[i]], file = fname, quote = TRUE, row.names = FALSE
    )
}

## Write predicted sets
for (i in seq_along(predictions_tips)) {
    fname <- paste0(
        gsub(' ', '_', phys_name), '_', rank_arg, '_phytools-ltp_predicted_Fold', i, '.csv'
    )
    write.csv(
        x = predictions_tips[[i]], file = fname, quote = TRUE, row.names = FALSE
    )
}

msg <- paste0('Script finihsed at ', Sys.time())
log_print(msg, blank_after = FALSE)
log_close()

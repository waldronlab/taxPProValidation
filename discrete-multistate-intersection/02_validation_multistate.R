args <- commandArgs(trailingOnly = TRUE)
suppressMessages({
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(mltools)
    library(logr)
})

phys_name <- args[[1]]
method <- args[[2]]

log_open(paste0(phys_name, '_', method, '_mcc'), logdir = TRUE, show_notes = TRUE)

dir <- file.path('.')

pattern <- paste0('^', phys_name, '_(all|genus|species|strain)_', method, '.*csv')
fnames <- list.files(path = dir, pattern = pattern, full.names = TRUE)
if (!length(fnames)) {
    msg <- paste0("Not enough data for prediction: ", phys_name)
    log_print(msg, blank_after = TRUE)
    quit(save = "no")
}

l <- map(fnames, read.csv)
names(l) <- sub("^(.*/)*(.*)\\.csv$", "\\2", fnames)
regex <- paste0('^', phys_name, "_(all|genus|species|strain)(.*)$")
new_phys_name <- gsub('_', ' ', phys_name)
names(l) <- paste0(new_phys_name, '_', sub(regex, '\\1\\2', names(l)))

test_folds <- l[grep('test', names(l))]
predicted_folds <- l[grep('predicted', names(l))]

myFun <- function(x, y) {
    x <- x |> 
        select(NCBI_ID, Attribute, Score) |> 
        complete(NCBI_ID, Attribute, fill = list(Score = 0)) |> 
        arrange(NCBI_ID, Attribute) |> 
        rename(tScore = Score) |> 
        mutate(tPN = ifelse(tScore > 0.8, 1, 0))
    y <- y |> 
        select(NCBI_ID, Attribute, Score) |> 
        group_by(NCBI_ID, Attribute) |> 
        slice_head(n = 1) |> 
        ungroup() |> 
        arrange(NCBI_ID, Attribute) |> 
        filter(NCBI_ID %in% unique(x$NCBI_ID)) |> 
        rename(pScore = Score) |> 
        mutate(pPN = ifelse(pScore >= 0.5, 1, 0))
    output <- left_join(x, y, by = c('NCBI_ID', 'Attribute')) |> 
        mutate(
            TF = case_when(
                tPN == 1 & pPN == 1 ~ 'TP',
                tPN == 1 & pPN == 0 ~ 'FN',
                tPN == 0 & pPN == 1 ~ 'FP',
                tPN == 0 & pPN == 0 ~ 'TN'
                
            )
        ) |> 
        {\(y) split(y, y$Attribute)}()
    return(output)
}

dat <- map2(test_folds, predicted_folds, myFun)
dat <- list_flatten(dat)

mcc <- map_dbl(dat, ~ mcc(preds = .x$pPN, actuals = .x$tPN))
names(mcc) <- sub('(test|predicted)_', '', names(mcc))

mcc_df <- data.frame(dat_name = names(mcc), mcc = unname(mcc)) |> 
    separate(
        col = 'dat_name', 
        into = c('physiology', 'rank', 'method', 'fold', 'attribute'),
        sep = '_'
    ) |> 
    group_by(
        physiology, rank, method, attribute
    ) |> 
    summarise(
        mean = mean(mcc, na.rm = TRUE),
        sd = sd(mcc, na.rm = TRUE)
    ) |> 
    ungroup()

fname <- paste0(phys_name, '_', method, '_mcc.tsv')
write.table(
    x = mcc_df, file = fname, sep = '\t', row.names = FALSE, quote = FALSE
)

log_close()

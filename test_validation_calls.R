## Run from the repo root: Rscript test_validation_calls.R
source("validation_calls.R")
stopifnot(
    identical(callPredicted(c(0.49, 0.5, 0.51, NA)), c(FALSE, TRUE, TRUE, FALSE)),
    identical(callTruth(c(0.8, 0.81, NA)), c(FALSE, TRUE, FALSE))
)
tie <- data.frame(NCBI_ID = c("g__1", "g__1", "g__2", "g__2"),
    Attribute = c("motility--FALSE", "motility--TRUE", "motility--TRUE", "motility--FALSE"),
    Score = c(0.5, 0.5, 0.3, 0.7), stringsAsFactors = FALSE)
out <- pickBinaryCall(tie)
stopifnot(nrow(out) == 2,
    out$Attribute[out$NCBI_ID == "g__1"] == "motility--TRUE",   # tie goes to TRUE
    out$Attribute[out$NCBI_ID == "g__2"] == "motility--FALSE")
## row order must not change the tie result
out2 <- pickBinaryCall(tie[c(2, 1, 4, 3), ])
stopifnot(identical(out2$Attribute, out$Attribute))
cat("ok\n")

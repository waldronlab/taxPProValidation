## Run from the repo root: Rscript test_attribute_scores.R
source("attribute_scores.R")
near <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-9))

## Unweighted Cohen's kappa by hand: po = 0.7, pe = 0.5*0.6 + 0.5*0.4 = 0.5
pred  <- c(rep("a", 5), rep("b", 5))
truth <- c("a", "a", "a", "a", "b", "a", "a", "b", "b", "b")
stopifnot(near(weighted_kappa(pred, truth, rep(1, 10), c("a", "b")), (0.7 - 0.5) / 0.5))
## equal weights of any size give the same kappa
stopifnot(near(weighted_kappa(pred, truth, rep(3, 10), c("a", "b")),
               weighted_kappa(pred, truth, rep(1, 10), c("a", "b"))))
## doubling the weight of a tip equals duplicating it
w <- c(2, rep(1, 9))
stopifnot(near(weighted_kappa(pred, truth, w, c("a", "b")),
               weighted_kappa(c(pred, pred[1]), c(truth, truth[1]), rep(1, 11), c("a", "b"))))
## two ordered states: weighted equals unweighted
stopifnot(near(weighted_kappa(pred, truth, rep(1, 10), c("a", "b"), positions = 1:2),
               weighted_kappa(pred, truth, rep(1, 10), c("a", "b"))))
## three ordered states: an adjacent error costs half an exact one.
## All predicted "m" while truth is lo/hi: unweighted po = 0, weighted po = 0.5
s <- c("lo", "m", "hi")
k_w <- weighted_kappa(c("m", "m", "lo", "hi"), c("lo", "hi", "lo", "hi"), rep(1, 4), s, 1:3)
k_u <- weighted_kappa(c("m", "m", "lo", "hi"), c("lo", "hi", "lo", "hi"), rep(1, 4), s)
stopifnot(k_w > k_u)
## a state that never occurs still sets the scale: with positions 1:4, an
## error between 1 and 2 costs a third
s4 <- c("1", "2", "3", "4")
k4 <- weighted_kappa(c("1", "2"), c("2", "2"), c(1, 1), s4, 1:4)
cm <- matrix(0, 4, 4); cm[1, 2] <- 0.5; cm[2, 2] <- 0.5
ag <- 1 - abs(outer(1:4, 1:4, "-")) / 3
ex <- sum(ag * outer(rowSums(cm), colSums(cm)))
stopifnot(near(k4, (sum(ag * cm) - ex) / (1 - ex)))
## chance agreement of 1 gives NA with a warning
stopifnot(is.na(suppressWarnings(weighted_kappa(c("a", "a"), c("a", "a"), c(1, 1), c("a", "b")))))

## Weights: identical distributions give weight 1; a band twice as common
## among targets gets twice the weight
stopifnot(near(nsti_weights(c(0.01, 0.03, 0.07), c(0.01, 0.03, 0.07)), c(1, 1, 1)))
w <- nsti_weights(c(0.01, 0.03), c(0.01, 0.01, 0.03))
stopifnot(near(w[1] / w[2], 2))

## R2 with equal weights is the squared correlation
x <- c(1, 2, 3, 4, 5); y <- c(1.1, 1.9, 3.2, 3.9, 5.3)
stopifnot(near(weighted_r2(y, x, rep(1, 5)), cor(x, y)^2))

## Lower bound ignores NA folds
stopifnot(near(lower_bound(c(0.8, 0.9, NA)), mean(c(0.8, 0.9)) - qt(0.975, 1) * sd(c(0.8, 0.9)) / sqrt(2)))
cat("ok\n")

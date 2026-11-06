# 0001. Assign cross-validation folds by species

- **Status:** Proposed
- **Date:** 2026-11-05
- **Deciders:** Levi Waldron, Claude Code

## Context and Problem Statement

Folds were assigned per tree tip. Strains of one species, and tips that share an NCBI id, landed in both training and test folds, so a taxon could be predicted from its own strains or from a duplicate of itself. That inflates MCC and R2, and the thresholds in waldronlab/bugphyzzExports ADR 0001 decide what is released, so inflated scores let unreliable predictions through.

## Decision

1. Strains and subspecies take their species' fold; taxa without a species take their own id as the key.
2. Numeric folds are assigned per NCBI id, not per tip.
3. Classes with fewer species than folds are dropped from that run, and a run with fewer than two classes quits with a message.
4. Guards check that no NCBI id is in both the training and test sets.

## Alternatives Considered

- **Keep per-tip folds and describe the leakage as a limitation.** Rejected: the scores set release thresholds.
- **Fold by genus.** Not adopted: too few genera for stable folds in most attributes, and it tests a harder problem than the export solves.

## Consequences

- Validation v1.0.5 reflects predictions for taxa whose species is unseen; some classes dropped below threshold (aerobic, spirillum, COGEM rating 3) and "tail" lost its row.
- Rare classes may disappear from the summary table, which drops their ASR rows at export.

References: #1, #3, #5.

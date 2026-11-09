# 0002. Call a predicted state positive at the score the export ships

- **Status:** Proposed
- **Date:** 2026-11-09
- **Deciders:** Levi Waldron

## Context and Problem Statement

The cross-validation turns each held-out prediction into a yes/no call before computing MCC. Since v1.0.2 it called a state positive only when its posterior score was above 0.8. The export (waldronlab/bugphyzzExports) keeps every imputed state with a score of 0.5 or more and drops "rarely" and "never" (`taxPPro::scores2Freq()`), and `importBugphyzz()` treats all kept rows as the taxon's state. So the MCC measured a stricter classifier than the one users receive.

A per-prediction calibration on the v1.0.5 and v1.0.6 folds showed the effect: classes with many mid-range scores (motility, extreme environment, animal pathogen) scored 0.05 to 0.12 lower under the 0.8 call than under the shipped rule, because true positives with scores between 0.5 and 0.8 counted as misses. Classes whose scores sit near 1 (Gram stain, spore formation) were unchanged.

## Decision

The 02 validation scripts call a predicted state positive when its score is 0.5 or more, the same rule the export applies. The held-out truth stays the curated state with a score above 0.8. The pass thresholds (MCC 0.7, R2 0.5; bugphyzzExports ADR 0001) do not change.

## Alternatives considered

- Keep the 0.8 call and have the export ship only scores above 0.8. Rejected: it would drop about a quarter of the imputed annotations, and the calibration shows precision for those rows is close to that of higher-scoring rows for the classes that pass.
- Use the score of each prediction as a per-row filter instead of a class-level pass. Rejected for now: above 0.8 the score carries almost no information about whether a prediction is right (aerobic precision is about 0.85 in every bin from 0.8 to 1), so a per-row cut-off removes correct and incorrect predictions alike.

## Consequences

- `validation_summary.tsv` describes the shipped data. Numbers are not comparable with v1.0.6 and earlier for discrete attributes; the README version log records the change.
- At rank `all` no class changes side of the 0.7 threshold, so the released data do not change because of this decision.
- The calibration script and results live in waldronlab/bugphyzz-review (`analysis/calibration/`) and should be re-run after each validation.

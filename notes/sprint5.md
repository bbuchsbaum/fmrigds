# Sprint 5 Plan — Meta-analysis Reducers & Evidence Combiners (of 8)

**Sprint:** 5 of 8  
**Duration:** 3 weeks  
**Prerequisites:** Sprints 1–4 complete (alignment, masking, optimizer, derivations, variance propagation).

---

## Objectives

1. Implement `reduce()` verb with fixed-effects (inverse-variance weighted) and random-effects (DerSimonian–Laird) meta-analysis across subjects.
2. Provide evidence combiners (Stouffer, Fisher) at the reduction level with weighting options.
3. Ensure df aggregation, weighting, and output metadata are correct.

---

## Scope & Tasks

### 1. `reduce()` Verb & Execution

**Files:** `R/verb-reduce.R` (new), `R/reduce-exec.R` (new), `R/compute.R`

- Define `reduce()` API: `method = c("fixed", "random", "stouffer", "fisher")`, `weights = c("1/var", "n_eff", "equal", "custom")`, optional grouping (aggregate by contrast or custom group).
- Execution logic:
  * Fixed-effects: weighted mean and variance using 1/var (with weight selection).
  * Random-effects: add between-subject variance via DerSimonian–Laird; allow fallback to fixed if tau² ≤ 0.
  * Evidence combiners: reuse Sprint 3 Stouffer/Fisher but at reduction scope.
- Handle missing subjects (NA) gracefully; allow `mask_policy` application before reduce.

### 2. Weights & DF Handling

- Implement weight calculators (`weights_fixed`, `weights_random`, `weights_custom`).
- Satterthwaite df aggregation for reducers; support custom df when provided.
- Ensure metadata records weight scheme, tau² estimates, degrees of freedom per output.

### 3. Integration with Compute Pipeline

- Extend `.apply_plan_nodes()` to process reduce nodes after map.
- Update output GDS object: single subject slice (subject dimension = 1) or drop as vector depending on design; set metadata (e.g., `metadata$reducer` details).

### 4. Documentation & Tests

- Roxygen docs for new verb/helpers; update NEWS.
- Tests (`tests/testthat/test-reduce.R`):
  * Fixed-effects aggregator vs analytical results.
  * Random-effects aggregator vs known tau² example.
  * Evidence-only pipeline (z-only) using Stouffer/Fisher.
  * Integration: gds -> derive -> map -> reduce -> compute end-to-end.

---

## Deliverables

- ✅ `reduce()` verb (fixed, random, stouffer, fisher) with weighting options.
- ✅ Execution logic integrated into `compute()`.
- ✅ Weight, tau², df calculations with metadata annotations.
- ✅ Comprehensive tests for meta-analysis scenarios.

---

## Deferred / Out-of-Scope

- Advanced random-effects estimators (REML) — Sprint 6/7.
- Bayesian meta-analysis (future work).
- Optimization of combiners via sparse operations — Sprint 8.

---

## Acceptance Criteria

1. Fixed-effects reductions match analytical inverse-variance weighted results.
2. Random-effects reductions compute tau², adjust weights, and report final variance/df.
3. Evidence-only reducers (Stouffer/Fisher) operate when beta/var absent.
4. Integration tests show alignment/mask/map → reduce pipeline.
5. Metadata captures method, weights, tau², df.
6. CI passes with new tests and documentation updates.

---

## Looking Ahead

- Sprint 6: HDF5 backend (gds-h5), plan serialization to disk, provenance append.
- Sprint 7: fmristore adapter integration, exporters (`write_out()`), factorial utilities.
- Sprint 8: Performance tuning, pkgdown, release readiness.

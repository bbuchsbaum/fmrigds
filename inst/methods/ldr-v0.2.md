# Local Displacement Rescue v0.2 whole-search contract

## Scope and estimand

LDR v0.2 extends the exact, one-center v0.1 reference to a prespecified search
set. With `centers = NULL`, that set is every active voxel whose complete
physical patch lies inside the image and active mask. A restricted `centers`
argument changes the multiplicity family and must be recorded in provenance.

The procedure flags a center only when all of the following hold:

1. a cross-fitted positive signed feature is reproducible there;
2. reproducible intersubject variation lies in the spatial-derivative tangent
   space of that feature;
3. both claims survive maximum-statistic correction over every searched center
   and patch/shift configuration;
4. an exact discrete-shift refit improves on both the aligned-positive and
   aligned signed-amplitude competitors; and
5. exact-model activation prevalence exceeds a prespecified lower bound; and
6. the ordinary fixed-coordinate group test does not survive its own
   maximum-statistic correction.

The output therefore attributes a missed result to bounded spatial discordance.
It does not attribute the discordance to anatomical registration error.

## Cross-fitted tangent screen

For each training fold, a compact signed template `g` is estimated from the
paired split maps. Finite translations of `g` define a physical derivative
basis `G`. The held-out split maps are fitted separately to `[B, g, G]`, where
`B` contains the local intercept and available linear spatial trends.

Activation evidence is the directional difference between positive- and
negative-amplitude soft matched-filter log evidence for the training template,
integrated over the bounded shift family with a fixed training-fold amplitude
scale. It does not select the best held-out shift, and reversing a held-out
map reverses its directional contribution. If the held-out tangent
coefficients are `b_is`, displacement evidence uses the cross-split product

```
D_i = b_iA' b_iB.
```

Before forming this product, the fitted tangent pattern is centered among the
held-out subjects within each fold, enforcing the identifying constraint
`E(delta_i) = 0` relative to that fold's training template and preventing a
common template-location error from being counted as displacement. Fold
assignment is a stable hash of subject identity, not acquisition or storage
order.
Independent split noise then has zero expected cross-product. A reproducible small
translation has `b_is` approximately equal to `-a_i G delta_i` and therefore
positive expected cross-split tangent energy. Jointly fitting `g` and `G`
prevents ordinary amplitude variation from being counted as displacement.

Local noise correlation is estimated only from training-subject standardized
split differences. The estimate is shrunk toward identity, made positive
definite, and used in held-out generalized least-squares fits. This covariance
affects the screen; the exact v0.1 confirmation remains the supplied-diagonal
likelihood until a separately validated correlated-noise likelihood is added.

## Search and calibration

Activation and ordinary fixed-coordinate statistics are calibrated with paired
whole-subject sign flips. Every flip reruns template estimation, tangent
construction, covariance estimation, all centers, all configurations, and
cross-fitting. The maximum statistic from each rescan provides strong
search-family control under the symmetric one-sample null.

The aligned-effect displacement null is calibrated by a Rademacher wild
bootstrap of the centered held-out `D_i` contributions. One multiplier is
shared by all centers and configurations for a subject, preserving their
observed spatial dependence. Every realization contributes the maximum
studentized displacement statistic. This is the predictive-score bootstrap
route; sign flips of the original activation maps are not used as a
displacement null.

For each component, the corrected p-value is

```
(1 + number of null maxima at least as large as observed) / (B + 1).
```

The screen conjunction is the maximum of corrected activation and displacement
p-values. Scale choice is paid for because null maxima include all requested
configurations. The ordinary MNI failure condition is evaluated with its own
corrected p-value and is not folded into the rescue p-value.

## Exact confirmation and fail-closed behavior

Only screen-significant centers are passed to the exact discrete-shift model.
Exact confirmation requires positive activation evidence, positive
improvements over both aligned competitors, positive MNI loss, and activation
prevalence at least `min_prevalence`. It can remove a screen discovery but
cannot create one, so it cannot increase the screen familywise error rate.

The method fails closed for incomplete patches, missing active-mask voxels,
insufficient subjects, rank-deficient feature/tangent designs, non-finite
precision inputs, and searches with no eligible centers. Non-eligible active
voxels receive missing inferential values and a zero LDR flag.
The tangent displacement scale is missing when too few subjects have stable
positive-amplitude linearized shift estimates or when cross-split shift
variance is nonpositive; the exact refit remains authoritative in that regime.

## Outputs

The principal `ldr_flag` equals `-log10(ldr_fwer_p)` only for exact-confirmed
screen discoveries whose ordinary MNI test fails; it is zero otherwise.
Companion maps report componentwise corrected p-values, tangent statistics,
MNI loss, displacement scale, and the fraction of reproducible feature energy
lying in the tangent directions. Full search configuration, resample maxima,
fold assignment, candidate receipts, and any exact-refit truncation are stored
as metadata.

## Remaining boundary

The v0.2 procedure is a translation-only voxel-space method. It does not cross
sulcal banks safely, estimate local diffeomorphisms, accept single-map inputs,
or identify the biological or preprocessing cause of displacement. Two nearby
functional areas that produce indistinguishable signed patches in different
subjects are observationally equivalent to a displaced common feature in these
inputs and cannot be separated without anatomy, surfaces, labels, or another
source of information.

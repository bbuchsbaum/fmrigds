# Local Displacement Rescue v0.1 contract

## Estimand

LDR v0.1 asks whether a prespecified MNI voxel anchors a compact signed local
feature that is reproducible across subjects and paired data splits, and whether
bounded subject-specific translations are needed to predict held-out subjects.
It estimates spatial-discordance rescue, not anatomical registration error.

## Reference models

For subject `i`, split `s`, and the patch anchored at `v`, the fitted family is

```
y_isv = B_v gamma_is + z_i a_i S_delta_i g_v + epsilon_isv.
```

The template has unit norm and a nonnegative anchor voxel; values outside the
anchor remain signed. The paired splits share `z_i`, `a_i`, and `delta_i`.
The reference implementation fits four models:

- `M0`: no local feature;
- `MA`: an aligned positive-amplitude feature;
- `MH`: an aligned signed-amplitude heterogeneity feature;
- `MJ`: a positive-amplitude feature with a centered discrete shift prior.

Displacement evidence is the smaller held-out improvement of `MJ` over `MA`
and `MH`. The `MH` comparison is necessary because a stable polarity mixture
can otherwise be represented by translations of an asymmetric signed template.

## Current output

`local_displacement_rescue()` is an exact one-center, one-contrast reference
implementation. It returns an optional ROI intersection-union p-value, a
model-based counterfactual change in the ordinary one-sample t statistic, and
the selected shift-prior RMS scale. Activation is calibrated with paired
whole-subject sign flips. Displacement and MNI loss are calibrated with
parametric bootstraps from both aligned competitors. Every resample reruns
cross-fitting and displacement-scale selection.

The current implementation uses supplied diagonal `var` or `se` precision,
removes a constant and available linear patch trends, integrates over shifts,
selects displacement scale only in training folds, and fails closed at image or
mask boundaries.

## Explicitly deferred

- spatial covariance estimated from split differences;
- multiple-comparison correction;
- tangent-space screening;
- whole-brain execution;
- single-map inputs and transformations beyond translation.

No corrected LDR flag may be exposed until all three component claims are
calibrated and the full adaptive search is included in each null realization.

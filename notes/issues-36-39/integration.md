# Issues 36-39: upstream integration and push preparation

Rebased the local fixes onto origin/main at
`e9acb89ff5796612fd17d01665da509cae1e81a8` (60 commits newer than the original
checkout). Final changed source/test/document hashes are in
`rebased-source-sha256.txt`; the earlier manifest and logs retain historical
pre-rebase evidence.

Resolved overlapping reducer changes by retaining upstream's corrected DL
meta-regression denominator, validity/rank checks, normal-inference outputs,
and OpenMP loading and fork controls. Applied the new no-floor default to the
new R regression fallback and the fixed/random examination diagnostic kernels.
Updated the upstream variance-floor tests to request the floor explicitly and
assert its warning, preserving their numeric and invalid-variance checks.
Updated NEWS and the new variance-floor documentation section.

Built native code from scratch in a temporary source copy, omitting the
checkout's tracked binaries. Full headless testthat run on the integrated code:

- 3855 assertions passed, including all 264 assertions in test-issues-36-39.R.
- 0 assertion failures and 0 test warnings.
- 13 test errors: 12 in test-frame-fit.R and 1 in test-frame-gds-projection.R,
  all because the installed multidesign does not export design_spec.
- 2 deliberate upstream skips: OMP fork test on macOS, and installed-package
  loader test during a pkgload development run.

The same frame-fit and frame-projection errors were reproduced against a fresh
native build of untouched upstream e9acb89. These errors are not new regressions;
the full integrated suite is not claimed as passing. Raw logs and exit-status
sidecars are retained here. Arrow sandbox CPU-info diagnostics and the testthat
R-build-version startup warning remain environment messages.

Source/documentation git diff --check passes (raw logs are retained byte-for-byte
and contain console whitespace). No unrelated documentation regeneration or
tracked native binaries are included. Full R CMD check remains unrun because
15 GiB free is below the inherited 20 GiB threshold. Hosted CI and merging are
not claimed. The unrelated UPSTREAM-ISSUES-from-fmrireg.md remains untracked.

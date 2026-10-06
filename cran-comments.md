## Submission notes

* First submission.
* `fmridataset`, `multidesign`, `neurotabs`, `neurothresh` and `albersdown`
  are optional (Suggests) and are available from
  <https://bbuchsbaum.r-universe.dev> (declared in `Additional_repositories`).
  All uses are guarded with `requireNamespace()` and the corresponding tests
  are skipped when they are absent.
* Before uploading, remove the `Remotes:` field from DESCRIPTION (it is used
  only by GitHub Actions dependency resolution).
* OpenMP kernels default to at most two threads; set
  `options(fmrigds.threads = n)` to change this.

## R CMD check results

Local check (R 4.3.3, Ubuntu 24.04, `_R_CHECK_LIMIT_CORES_=true`):
0 errors | 0 warnings | 2 notes

* New submission / Suggests not available for checking: the r-universe
  packages listed above.
* Installed package size: about 4 MB of `libs` (RcppArmadillo kernels) and
  about 5 MB of `doc`, because each HTML vignette embeds the bundled web fonts.

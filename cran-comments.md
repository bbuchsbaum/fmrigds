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

0 errors | 0 warnings | 1 note

* New submission.

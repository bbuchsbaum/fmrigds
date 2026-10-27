# Run from the repository root against a freshly built temporary source copy:
# Rscript --vanilla scripts/repro-issues-36-39.R /path/to/source-copy
args <- commandArgs(trailingOnly = TRUE)
pkgload::load_all(if (length(args)) args[[1]] else ".", quiet = TRUE)
g <- new_gds(
  list(beta = array(c(2, NA_real_), c(2, 1, 1)),
       var = array(1, c(2, 1, 1))),
  space_sample_labels(c("kept", "unused")), "s1", "c1")
print(c(issue36_beta = assay(compute(map_to(
  g, space_sample_labels("target"), matrix(c(1, 0), 1))), "beta")[1]))
a <- diag(4)
b <- a
b[1, 4] <- 10
print(tryCatch(assert_compatible_spaces(
  space_voxel(c(2, 2, 2), a), space_voxel(c(2, 2, 2), b)),
  error = conditionMessage))
g <- new_gds(
  list(beta = array(c(1, 2), c(1, 1, 2)), var = array(1, c(1, 1, 2))),
  space_sample_labels("f1"), "s1", c("a", "b"))
print(tryCatch(contrasts(harmonise_contrasts(g, c(a = "b"))),
               error = conditionMessage))
fit <- function(scale = 1) {
  g <- new_gds(
    list(beta = array(c(1, 2, 3) * scale, c(1, 3, 1)),
         var = array(scale^2, c(1, 3, 1))),
    space_sample_labels("f1"), c("s1", "s2", "s3"), "c1")
  out <- compute(reduce(g, method = "fixed"))
  c(z = assay(out, "z_g")[1], p = assay(out, "p_g")[1])
}
print(rbind(original = fit(), rescaled = fit(1e-8)))
sessionInfo()

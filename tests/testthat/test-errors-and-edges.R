test_that("posthoc errors on unknown method", {
  set.seed(1)
  df <- data.frame(
    sample = rep(paste0("R",1:3), each = 2),
    subject = rep(c("s1","s2"), times = 3),
    contrast = "c1",
    p = runif(6)
  )
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE)
  plan <- gds(tmp, effect_cols = list(p = "p")) |> posthoc("not:a:method")
  expect_error(compute(plan), "Unknown post-hoc method")
})

test_that("use_weight errors when assay missing and attach_weight dims must match", {
  beta <- array(0, dim = c(2,2,1)); var <- array(1, dim = c(2,2,1))
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("a","b")), c("s1","s2"), "c1")
  bad_w <- array(1, dim = c(3,2,1))
  expect_error(attach_weight(g, "w", bad_w), "dims must match")
  expect_error(use_weight(g, "w_missing"), "Weight assay not found")
})

test_that("space_subset warns for dense voxel spaces", {
  sp <- space_voxel(dim = c(2,2,1), affine = diag(4), storage = "dense")
  expect_warning(space_subset(sp, 1:2), "cannot be reduced")
})

test_that("register_alignment refuses duplicate without overwrite", {
  from <- space_parcels(c("a","b")); to <- from
  ops <- list(s1 = diag(2), s2 = diag(2))
  fam <- make_linear_family("dup", from, to, ops)
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "a,s1,c1,1,1","b,s1,c1,1,1","a,s2,c1,1,1","b,s2,c1,1,1"
  ), tmp)
  plan <- gds(tmp)
  plan <- register_alignment(plan, fam)
  expect_error(register_alignment(plan, fam), "already registered")
})


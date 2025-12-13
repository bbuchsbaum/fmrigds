test_that("apply_align maps each subject into group space", {
  beta <- array(c(1, 2, 3, 4), c(2, 2, 1))
  var <- array(c(0.25, 0.5, 0.5, 1), c(2, 2, 1))
  arrays <- list(beta = beta, var = var, se = sqrt(var), t = beta / sqrt(var))

  S1 <- diag(2)
  S2 <- matrix(c(0, 1, 1, 0), nrow = 2)  # swap
  family <- MapFamily(
    name = "swap",
    from_space = space_parcels(c("roi1", "roi2")),
    to_space = space_parcels(c("grp1", "grp2")),
    type = "linear",
    by_subject = list(s1 = S1, s2 = S2)
  )

  res <- apply_align(family, arrays, c("s1", "s2"), family$from)
  beta_out <- res$arrays$beta

  expect_equal(beta_out[, 1, 1], c(1, 2))
  expect_equal(beta_out[, 2, 1], c(4, 3))
  expect_s3_class(res$space, "space_parcels")
})

test_that("map families register and persist through h5 write/read", {
  skip_if_not_installed("hdf5r")

  beta <- array(seq_len(8), c(4, 2, 1))
  var <- array(1, dim(beta))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )

  family <- MapFamily(
    name = "identity_parcel",
    from_space = g$space,
    to_space = g$space,
    type = "linear",
    by_subject = list(
      s1 = diag(4),
      s2 = diag(4)
    )
  )

  g <- register_map(g, family)
  expect_equal(list_map_families(g), "identity_parcel")

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  write_gds_h5(g, tmp)

  plan <- gds(tmp)
  expect_equal(list_map_families(plan), "identity_parcel")
  fam_loaded <- get_map_family(plan, "identity_parcel")
  expect_s3_class(fam_loaded, "gds_map_family")

  aligned <- compute(align(plan, "identity_parcel"))
  expect_equal(dim(assay(aligned, "beta")), c(4, 2, 1))
})

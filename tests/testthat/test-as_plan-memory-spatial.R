test_that("as_plan.gds enables verbs on realized GDS and se->var derivation", {
  set.seed(42)
  nS <- 4; nB <- 3; nK <- 1
  beta <- array(rnorm(nB * nS * nK), dim = c(nB, nS, nK))
  se   <- array(runif(nB * nS * nK, min = 0.1, max = 0.3), dim = c(nB, nS, nK))
  g <- as_gds(list(beta = beta, se = se),
              space = space_sample_labels(paste0("ROI_", seq_len(nB))),
              subjects = paste0("s", seq_len(nS)),
              contrasts = paste0("c", seq_len(nK)))

  # Use realized GDS directly with reduce(); var must be auto-derived from se
  gout <- reduce(g, method = "meta:fe") |> compute()

  # Expect group-level results mapped to standard assay names
  expect_true(all(c("beta", "var", "z", "p") %in% names(assays(gout))))
  # se_g should be present and consistent with var
  expect_true("se_g" %in% names(assays(gout)))
  vg <- as.numeric(assay(gout, "var")[, 1, 1])
  seg <- as.numeric(assay(gout, "se_g")[, 1, 1])
  expect_equal(vg, seg^2, tolerance = 1e-8)
})

test_that("posthoc fdr:spatial provides q and broadcasts within groups", {
  skip_if(!("fdr:spatial" %in% list_posthoc()))
  set.seed(24)
  nS <- 2; nB <- 6; nK <- 1
  p <- array(runif(nB * nS * nK), dim = c(nB, nS, nK))
  grp <- c(1, 1, 2, 2, 3, 3)
  g <- as_gds(list(p = p),
              space = space_sample_labels(paste0("ROI_", seq_len(nB))),
              subjects = paste0("s", seq_len(nS)),
              contrasts = paste0("c", seq_len(nK)))
  gout <- posthoc(g, method = "fdr:spatial", options = list(group = grp)) |> compute()
  expect_true("q" %in% names(assays(gout)))
  q <- assay(gout, "q")[, 1, 1]
  # Within-group q's should be equal (broadcast from group-level adjustment)
  expect_equal(q[1], q[2])
  expect_equal(q[3], q[4])
  expect_equal(q[5], q[6])
})


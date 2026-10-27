test_that("#36 zero map coefficients isolate missing effects and variances", {
  M <- rbind(c(2, 0), c(0, -1), c(0.5, 0.5), c(0, 0))
  for (missing in c(NA_real_, NaN, Inf)) {
    for (field in c("beta", "var")) {
      arrays <- list(beta = array(c(2, 3), c(2, 2, 2)),
                     var = array(c(4, 9), c(2, 2, 2)))
      arrays[[field]][2, , ] <- missing
      g <- new_gds(arrays, space_sample_labels(c("a", "b")),
                   c("s1", "s2"), c("c1", "c2"))
      for (map in list(M, Matrix::Matrix(M, sparse = TRUE))) {
        out <- compute(map_to(g, space_sample_labels(letters[1:4]), map))
        expect_equal(as.numeric(assay(out, "beta")[1, , ]), rep(4, 4))
        expect_equal(as.numeric(assay(out, "var")[1, , ]), rep(16, 4))
        expect_true(all(!is.finite(assay(out, field)[2:3, , ])))
        expect_equal(as.numeric(assay(out, field)[4, , ]), rep(0, 4))
      }
    }
  }
})

test_that("#36 contributing missing variances propagate through t and df", {
  arrays <- list(beta = array(c(2, NA), c(2, 1, 1)),
                 var = array(c(4, NA), c(2, 1, 1)),
                 t = array(c(1, NA), c(2, 1, 1)),
                 df = array(c(10, NA), c(2, 1, 1)))
  g <- new_gds(arrays, space_sample_labels(c("a", "b")), "s1", "c1")
  out <- compute(map_to(g, g$space, diag(2)))
  expect_equal(as.numeric(assay(out, "t")), c(1, NA))
  expect_equal(as.numeric(assay(out, "df")), c(10, NA))
  family <- MapFamily("identity", g$space, g$space, by_subject = list(s1 = diag(2)))
  aligned <- compute(align(register_map(g, family), "identity"))
  expect_equal(assay(aligned, "beta"), assay(out, "beta"))
  expect_equal(assay(aligned, "df"), assay(out, "df"))
})

test_that("#37 voxel compatibility includes the affine by default", {
  a <- diag(4)
  translation <- a; translation[1, 4] <- 10
  rotation <- a; rotation[1:2, 1:2] <- matrix(c(0, 1, -1, 0), 2)
  spacing <- a; spacing[1, 1] <- 2
  for (storage in c("dense", "packed")) {
    s1 <- space_voxel(c(2, 2, 2), a, mask_idx = 1:4, storage = storage)
    expect_true(assert_compatible_spaces(s1, s1))
    for (affine in list(translation, rotation, spacing)) {
      s2 <- space_voxel(c(2, 2, 2), affine, mask_idx = 1:4, storage = storage)
      expect_error(assert_compatible_spaces(s1, s2), "affine")
      expect_true(assert_compatible_spaces(s1, s2, fields = c("type", "dim")))
    }
  }
  expect_true(assert_compatible_spaces(space_parcels("a"), space_parcels("a")))
})

test_that("#38 contrast renaming rejects collisions and preserves metadata", {
  beta <- array(c(1, 2), c(1, 1, 2), list("f1", "s1", c("a", "b")))
  g <- new_gds(list(beta = beta, var = beta^2),
               space_sample_labels("f1"), "s1", c("a", "b"))
  g <- with_contrast_data(g, data.frame(value = c(20, 10), row.names = c("b", "a")))
  expect_error(harmonise_contrasts(g, c(a = "b")), "unique")
  expect_error(harmonise_contrasts(g, c(a = "x", b = "x")), "unique")
  expect_error(harmonise_contrasts(g, c(a = "x", a = "y")), "unique")
  expect_error(harmonise_contrasts(g, c(a = NA_character_)), "nonempty")
  expect_error(harmonise_contrasts(g, c(a = "")), "nonempty")
  expect_equal(contrasts(g), c("a", "b"))
  for (mapping in list(c(a = "A"), c(a = "b", b = "a"))) {
    out <- harmonise_contrasts(g, mapping)
    expect_equal(as.numeric(assay(out, "beta")), c(1, 2))
    expect_equal(dimnames(assay(out, "beta"))[[3]], contrasts(out))
    expect_equal(dimnames(assay(out, "var"))[[3]], contrasts(out))
    expect_equal(rownames(contrast_data(out)), contrasts(out))
    expect_equal(contrast_data(out)$value, c(10, 20))
    one <- compute(subset(as_plan(out), contrast = contrasts(out)[2]))
    expect_equal(as.numeric(assay(one, "beta")), 2)
    expect_equal(contrast_data(one)$value, 20)
  }
})

# Exercise the actual R fallback body without changing the package namespace.
meta_r_fallback <- function(fun) {
  env <- new.env(parent = environment(fun))
  env$exists <- function(...) FALSE
  environment(fun) <- env
  fun
}

test_that("#39 fixed and random meta-analysis preserve inference across units", {
  beta <- matrix(c(1, 2, 3, -2, 1, 7), 3)
  var <- matrix(c(1, 1, 1, 1, 2, 4), 3)
  kernels <- list(core_meta_fe_kernel, core_meta_re_dl_kernel,
                  meta_r_fallback(core_meta_fe_kernel),
                  meta_r_fallback(core_meta_re_dl_kernel))
  for (kernel in kernels) {
    base <- kernel(beta, var)
    # Equal-variance analytic case: mean = 2, variance of mean = 1/3.
    expect_equal(base$beta_g[1], 2)
    expect_equal(base$var_g[1], 1 / 3)
    expect_equal(base$z_g[1], 2 * sqrt(3))
    expect_equal(base$p_g[1], 2 * pnorm(-2 * sqrt(3)))
    expect_equal(base$Q[2], 117 / 7)
    if ("tau2" %in% names(base)) {
      expect_equal(base$tau2[2], 103 / 7)
    } else {
      expect_equal(base$beta_g[2], 1 / 7)
      expect_equal(base$var_g[2], 4 / 7)
    }
    for (scale in c(1e-8, 1e8)) {
      out <- kernel(beta * scale, var * scale^2)
      for (name in c("z_g", "p_g", "Q", "I2", "n_eff")) {
        expect_equal(out[[name]], base[[name]], tolerance = 1e-10, info = name)
      }
      for (name in c("beta_g", "se_g")) {
        expect_equal(out[[name]] / scale, base[[name]], tolerance = 1e-10)
      }
      for (name in intersect(c("var_g", "tau2"), names(base))) {
        expect_equal(out[[name]] / scale^2, base[[name]], tolerance = 1e-10)
      }
    }
    # Invalid observations must not be made valid by a zero or positive floor.
    y <- rbind(beta, 100, 100, 100, NA, Inf)
    v <- rbind(var, 0, -1, Inf, 1, 1)
    expect_equal(kernel(y, v), base)
    expect_warning(kernel(beta * 1e-8, var * 1e-16, opts = list(eps = 1e-12)),
                   "Explicit eps variance floor applied")
    expect_error(kernel(beta, var, opts = list(eps = -1)), "nonnegative")
  }
})

test_that("#39 public fixed reduction retains the reported analytic z", {
  for (scale in c(1, 1e-8)) {
    g <- new_gds(list(beta = array((1:3) * scale, c(1, 3, 1)),
                      var = array(scale^2, c(1, 3, 1))),
                 space_sample_labels("f1"), paste0("s", 1:3), "c1")
    out <- compute(reduce(g, method = "fixed"))
    expect_equal(as.numeric(assay(out, "z_g")), 2 * sqrt(3))
  }
})

test_that("#39 meta-regression matches intercept-only meta-analysis across units", {
  beta <- matrix(c(-2, 1, 7, 5), 4)
  var <- matrix(c(1, 2, 4, 3), 4)
  X <- matrix(1, 4, 1)
  for (method in c("fe", "re")) {
    cpp <- if (method == "fe") meta_fe_reg_cpp else meta_re_reg_dl_cpp
    r <- get_reducer(paste0("meta:", method, "_reg"))$fun
    regressions <- list(cpp, function(beta, var, X) r(beta, var, X, opts = list()))
    if (method == "re") regressions <- c(regressions, list(.meta_re_reg_dl_r))
    kernel <- if (method == "fe") core_meta_fe_kernel else core_meta_re_dl_kernel
    base <- kernel(beta, var)
    for (scale in c(1, 1e-8, 1e8)) {
      for (reg in regressions) {
        out <- reg(beta * scale, var * scale^2, X)
        expect_equal(as.numeric(out$coef) / scale, base$beta_g, tolerance = 1e-10)
        expect_equal(as.numeric(out$se_coef) / scale, base$se_g, tolerance = 1e-10)
        expect_equal(out$Q, base$Q, tolerance = 1e-10)
        if (method == "re") expect_equal(out$tau2 / scale^2, base$tau2, tolerance = 1e-10)
        invalid <- reg(rbind(beta * scale, 100, 100, 100),
                       rbind(var * scale^2, 0, -1, Inf), rbind(X, 1, 1, 1))
        expect_equal(invalid, out)
      }
    }
  }
})

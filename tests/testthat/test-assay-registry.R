test_that("register_assay() registers assays correctly", {
  # Register a custom assay
  result <- register_assay("custom_stat", role = "z", units = "zscore")

  expect_equal(result$name, "custom_stat")
  expect_equal(result$role, "z")
  expect_equal(result$units, "zscore")
})

test_that("register_assay() validates role argument", {
  expect_error(
    register_assay("invalid", role = "not_a_role"),
    "'arg' should be one of"
  )
})

test_that("register_assay() handles variance_of parameter", {
  result <- register_assay("custom_var", role = "variance", variance_of = "custom_beta")
  expect_equal(result$variance_of, "custom_beta")
})

test_that("register_assay() handles derive_from parameter", {
  result <- register_assay("custom_t", role = "t", derive_from = c("beta", "var", "df"))
  expect_equal(result$derive_from, c("beta", "var", "df"))
})

test_that("assay_info() retrieves registered assays", {
  # Use a default assay
  info <- assay_info("beta")
  expect_equal(info$name, "beta")
  expect_equal(info$role, "location")
  expect_equal(info$units, "%BOLD")
})

test_that("assay_info() returns NULL for unknown assays", {
  info <- assay_info("nonexistent_assay")
  expect_null(info)
})

test_that("can_map_linear() returns TRUE for mappable assays", {
  expect_true(can_map_linear("beta"))   # location
  expect_true(can_map_linear("var"))    # variance
  expect_true(can_map_linear("se"))     # stdev
})

test_that("can_map_linear() returns FALSE for non-mappable assays", {
  expect_false(can_map_linear("t"))     # t-statistic
  expect_false(can_map_linear("z"))     # z-statistic
  expect_false(can_map_linear("p"))     # p-value
  expect_false(can_map_linear("F"))     # F-statistic
})

test_that("can_map_linear() returns FALSE for unknown assays", {
  expect_false(can_map_linear("unknown_assay"))
})

test_that("default assays are registered on package load", {
  # These should all exist
  expect_true(!is.null(assay_info("beta")))
  expect_true(!is.null(assay_info("var")))
  expect_true(!is.null(assay_info("se")))
  expect_true(!is.null(assay_info("t")))
  expect_true(!is.null(assay_info("z")))
  expect_true(!is.null(assay_info("F")))
  expect_true(!is.null(assay_info("df")))
  expect_true(!is.null(assay_info("df1")))
  expect_true(!is.null(assay_info("df2")))
  expect_true(!is.null(assay_info("n_eff")))
  expect_true(!is.null(assay_info("p")))
  expect_true(!is.null(assay_info("chi2")))
  expect_true(!is.null(assay_info("logBF")))
})

test_that("assay roles are correct for default assays", {
  expect_equal(assay_info("beta")$role, "location")
  expect_equal(assay_info("var")$role, "variance")
  expect_equal(assay_info("se")$role, "stdev")
  expect_equal(assay_info("t")$role, "t")
  expect_equal(assay_info("z")$role, "z")
  expect_equal(assay_info("F")$role, "F")
  expect_equal(assay_info("df")$role, "df")
  expect_equal(assay_info("p")$role, "p")
  expect_equal(assay_info("chi2")$role, "chi2")
  expect_equal(assay_info("logBF")$role, "log_evidence")
})

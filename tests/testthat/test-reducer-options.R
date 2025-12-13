test_that("reducer options schema validates enums", {
  # Access internal helper directly
  vo <- fmrigds:::validate_reducer_options
  # Valid value passes
  expect_equal(vo(list(return_cov = c("none", "tri")), list(return_cov = "tri"))$return_cov, "tri")
  # Missing defaults to first allowed
  expect_equal(vo(list(return_cov = c("none", "tri")), list())$return_cov, "none")
  # Invalid value errors
  expect_error(vo(list(return_cov = c("none", "tri")), list(return_cov = "bad")), "Invalid option 'return_cov'")
})


test_that("subset_eager() and derive_eager() exist and are callable", {
  # These functions exist and wrap compute()
  # Testing them requires full adapter setup which is covered by other tests
  # Here we just verify they exist and have the right signature

  expect_true(is.function(subset_eager))
  expect_true(is.function(derive_eager))

  # Verify they call compute() by checking their body
  subset_eager_body <- paste(deparse(body(subset_eager)), collapse = " ")
  derive_eager_body <- paste(deparse(body(derive_eager)), collapse = " ")

  expect_true(grepl("compute", subset_eager_body, fixed = TRUE))
  expect_true(grepl("compute", derive_eager_body, fixed = TRUE))
})

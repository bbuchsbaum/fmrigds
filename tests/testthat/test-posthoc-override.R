test_that("register_posthoc overwrite guard and replacement work", {
  name <- "tmp:method"
  on.exit({ try(fmrigds::unregister_posthoc(name), silent = TRUE) }, add = TRUE)
  f1 <- function(arrays, opts) list(q = arrays$p)
  f2 <- function(arrays, opts) list(q = 1 - arrays$p)
  register_posthoc(name, f1, requires = c("p"), provides = c("q"), overwrite = FALSE)
  expect_error(register_posthoc(name, f2, requires = c("p"), provides = c("q"), overwrite = FALSE))
  register_posthoc(name, f2, requires = c("p"), provides = c("q"), overwrite = TRUE)
  ph <- get_posthoc(name)
  expect_true(identical(ph$fun, f2))
})


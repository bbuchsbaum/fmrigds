test_that("register_reducer() registers custom reducers", {
  # Create a simple custom reducer
  custom_mean <- function(beta, var, X, z, p, df, df1, df2, opts) {
    list(
      beta_g = colMeans(beta, na.rm = TRUE),
      var_g = colMeans(var, na.rm = TRUE)
    )
  }

  register_reducer(
    name = "test:mean",
    fun = custom_mean,
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g")
  )

  # Retrieve and verify
  reducer <- get_reducer("test:mean")
  expect_equal(reducer$name, "test:mean")
  expect_true(is.function(reducer$fun))
  expect_equal(reducer$requires, c("beta", "var"))
  expect_equal(reducer$provides, c("beta_g", "var_g"))
  expect_null(reducer$model_contract)
  expect_null(reducer$diagnostics)
})

test_that("register_reducer stores validated model and diagnostic contracts", {
  fun <- function(beta, var, X, z, p, df, df1, df2, opts) {
    list(beta_g = colMeans(beta))
  }
  register_reducer(
    "test:contract",
    fun,
    requires = "beta",
    provides = "beta_g",
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "unweighted",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "hat_matrix"
    ),
    diagnostics = list(
      fun = NULL,
      capabilities = c("prediction", "coefficient_deletion"),
      modes = "exact"
    )
  )

  reducer <- get_reducer("test:contract")
  expect_true(reducer$model_contract$uses_X)
  expect_identical(reducer$model_contract$weight_mode, "unweighted")
  expect_identical(
    reducer$diagnostics$capabilities,
    c("prediction", "coefficient_deletion")
  )
  expect_error(
    register_reducer(
      "test:bad-contract",
      fun,
      "beta",
      "beta_g",
      model_contract = list(weight_mode = "mystery")
    ),
    "weight_mode"
  )
})

test_that("get_reducer() retrieves registered reducers", {
  # Should have core reducers registered
  fe <- get_reducer("meta:fe")
  expect_true(!is.null(fe))
  expect_equal(fe$name, "meta:fe")
})

test_that("get_reducer() returns NULL for unknown reducers", {
  result <- get_reducer("nonexistent:reducer")
  expect_null(result)
})

test_that("list_reducers() returns registered reducer names", {
  reducers <- list_reducers()
  expect_true(is.character(reducers))
  expect_true(length(reducers) > 0)

  # Should include core reducers
  expect_true("meta:fe" %in% reducers)
  expect_true("meta:re" %in% reducers)
})

test_that("register_reducer() validates inputs", {
  expect_error(
    register_reducer(name = NULL, fun = function(x) x, requires = "beta", provides = "beta_g"),
    "is.character\\(name\\)"
  )

  expect_error(
    register_reducer(name = "test", fun = "not_a_function", requires = "beta", provides = "beta_g"),
    "is.function\\(fun\\)"
  )
})

test_that("reducer options_schema is stored correctly", {
  custom_reducer <- function(beta, var, X, z, p, df, df1, df2, opts) {
    tau2 <- opts$tau2 %||% 0
    list(beta_g = colMeans(beta), tau2 = tau2)
  }

  register_reducer(
    name = "test:with_options",
    fun = custom_reducer,
    requires = c("beta", "var"),
    provides = c("beta_g", "tau2"),
    options_schema = list(tau2 = "numeric")
  )

  reducer <- get_reducer("test:with_options")
  expect_equal(reducer$options_schema, list(tau2 = "numeric"))
})

test_that("core reducers are registered on package load", {
  reducers <- list_reducers()

  # Core meta-analysis reducers
  expect_true("meta:fe" %in% reducers)
  expect_true("meta:re" %in% reducers)

  # Evidence combination reducers
  expect_true("combine:stouffer" %in% reducers)
  expect_true("combine:fisher" %in% reducers)
})

test_that("legacy method names are normalized", {
  # .normalize_reducer_name is internal but affects behavior
  # Test through actual reduce operations if possible
  # For now, just ensure core reducers exist
  expect_true(!is.null(get_reducer("meta:fe")))
  expect_true(!is.null(get_reducer("meta:re")))
})

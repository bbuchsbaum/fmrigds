# Regression tests for engine review fixes (map/align assay propagation,
# reducer dispatch, weights, block/preview axes, plan serialization, subset
# on voxel spaces and realised objects, plan axis accessors, Lancaster df).


.review_stat_gds <- function(n = 6L, n_subj = 3L, df_value = 20, na_first = TRUE,
                             space = NULL, row_data = NULL) {
  set.seed(42)
  beta <- array(stats::rnorm(n * n_subj), dim = c(n, n_subj, 1L))
  var <- array(stats::runif(n * n_subj, 0.1, 0.5), dim = c(n, n_subj, 1L))
  df <- array(df_value, dim = c(n, n_subj, 1L))
  if (na_first) beta[1L, 1L, 1L] <- NA_real_
  t <- beta / sqrt(var)
  p <- 2 * stats::pt(-abs(t), df)
  z <- stats::qnorm(1 - p / 2) * sign(t)
  q <- array(0.5, dim = dim(beta))
  new_gds(
    assays = list(beta = beta, var = var, df = df, se = sqrt(var), t = t, z = z, p = p, q = q),
    space = space %||% space_sample_labels(paste0("r", seq_len(n))),
    subjects = paste0("s", seq_len(n_subj)),
    contrasts = "c1",
    row_data = row_data
  )
}

.expect_stats_consistent <- function(out) {
  a <- assays(out)
  expect_equal(a$t, a$beta / sqrt(a$var))
  expect_equal(a$se, sqrt(a$var))
  if (!is.null(a$df)) {
    expect_equal(a$p, 2 * stats::pt(-abs(a$t), a$df))
  } else {
    expect_equal(a$p, 2 * stats::pnorm(-abs(a$z)))
  }
  expect_equal(a$p, 2 * stats::pnorm(-abs(a$z)), tolerance = 1e-6)
}

# 1. map_to()/align() ------------------------------------------------------

test_that("mask() |> map_to() re-derives p/z/t and drops stale assays", {
  g <- .review_stat_gds()
  S <- matrix(c(0.5, 0.5, 0, 0, 0,
                0, 0, 1 / 3, 1 / 3, 1 / 3), nrow = 2, byrow = TRUE)
  target <- space_sample_labels(c("A", "B"))

  out <- g |> mask() |> map_to(target, map = S) |> compute()

  a <- assays(out)
  expect_false("q" %in% names(a))
  expect_true(all(vapply(a, function(x) dim(x)[1L], integer(1)) == 2L))
  expect_equal(dim(a$df), c(2L, 3L, 1L))
  .expect_stats_consistent(out)

  # FDR post-hoc consumes the updated p, not the stale source p.
  fdr <- g |> mask() |> map_to(target, map = S) |> posthoc("fdr:bh") |> compute()
  pa <- assay(fdr, "p")
  expect_equal(assay(fdr, "q")[, 1, 1], stats::p.adjust(pa[, 1, 1], method = "BH"))
  expect_equal(pa, a$p)
})

test_that("map_to() errors clearly when the map does not match the current sample axis", {
  g <- .review_stat_gds()
  S6 <- matrix(1 / 6, nrow = 1, ncol = 6)
  expect_error(
    g |> mask() |> map_to(space_sample_labels("A"), map = S6) |> compute(),
    "6 columns but the data currently has 5 samples"
  )
})

test_that("map_to() with df_rule = 'none' drops non-constant df", {
  g <- .review_stat_gds(na_first = FALSE)
  g$assays$df[2L, 1L, 1L] <- 5
  S <- matrix(c(1, 1, 0, 0, 0, 0), nrow = 1) / 2
  out <- g |>
    map_to(space_sample_labels("A"), map = S, uncertainty = UncertaintyRule("independent", df_rule = "none")) |>
    compute()
  a <- assays(out)
  expect_null(a$df)
  expect_equal(a$z, a$beta / sqrt(a$var))
  .expect_stats_consistent(out)
})

test_that("align() keeps df on the target sample axis and re-derives statistics", {
  g <- .review_stat_gds(na_first = FALSE)
  ops <- lapply(1:3, function(j) matrix(c(rep(1 / 3, 3), 0, 0, 0, 0, 0, 0, rep(1 / 3, 3)), nrow = 2, byrow = TRUE) * j)
  names(ops) <- subjects(g)
  fam_none <- MapFamily("fam_none", g$space, space_sample_labels(c("A", "B")),
                        by_subject = ops,
                        uncertainty = UncertaintyRule("independent", df_rule = "none"))
  out <- g |> align(fam_none) |> compute()
  a <- assays(out)
  expect_false("q" %in% names(a))
  expect_equal(dim(a$df), c(2L, 3L, 1L))
  expect_true(all(a$df == 20))
  .expect_stats_consistent(out)

  fam_sw <- MapFamily("fam_sw", g$space, space_sample_labels(c("A", "B")), by_subject = ops)
  out_sw <- g |> align(fam_sw) |> compute()
  expect_equal(dim(assay(out_sw, "df")), c(2L, 3L, 1L))
  .expect_stats_consistent(out_sw)
})

# 2./3. reduce() dispatch and weights ---------------------------------------

test_that("reduce() errors on unknown reducers instead of running Fisher", {
  g <- .review_stat_gds(na_first = FALSE)
  expect_error(reduce(g, "meta:RE"), "Unknown reducer 'meta:RE'.*meta:re")
  node <- op_reduce("not-a-reducer", "1/var", "contrast")
  expect_error(fmrigds:::apply_reduce(node, assays(g), "1/var", subjects(g)), "Unknown reducer")
})

test_that("reduce() refuses weight schemes the reducer would ignore", {
  g <- .review_stat_gds(na_first = FALSE)
  expect_error(reduce(g, "fixed", weights = "equal"), "not supported by reducer 'meta:fe'")
  expect_error(reduce(g, "fixed", weights = "n_eff"), "not supported")
  expect_error(reduce(g, "random", weights = "custom"), "not supported")
  w <- use_weight(attach_weight(g, "w", array(1, dim(assay(g, "beta")))), "w")
  expect_error(reduce(g, "fixed", weights = w$weights, options = w$options), "not supported")
  # Equal weighting is what unweighted / evidence reducers already do.
  expect_s3_class(reduce(g, "stouffer", weights = "equal"), "gds_plan")
  expect_s3_class(reduce(g, "ols:voxelwise", weights = "equal", formula = ~ 1), "gds_plan")
})

test_that("reducers declaring weight schemes receive them", {
  env <- fmrigds:::.gds_reducers
  on.exit(if (exists("test:review_weighted", envir = env)) rm("test:review_weighted", envir = env), add = TRUE)
  seen <- new.env()
  register_reducer(
    "test:review_weighted",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      seen$weights <- opts$weights
      list(beta_g = colMeans(beta), var_g = rep(1, ncol(beta)))
    },
    requires = c("beta"),
    provides = c("beta_g", "var_g"),
    options_schema = list(weights = c("1/var", "equal"))
  )
  g <- .review_stat_gds(na_first = FALSE)
  out <- reduce(g, "test:review_weighted", weights = "equal") |> compute()
  expect_equal(seen$weights, "equal")
  expect_equal(subjects(out), "meta")
  expect_error(reduce(g, "test:review_weighted", weights = "n_eff"), "not supported")
})

# 4. preview()/compute(block =) ---------------------------------------------

test_that("preview() subsets sample_labels space and row_data to the block", {
  rd <- data.frame(roi = paste0("roi", 1:10), stringsAsFactors = FALSE)
  g <- .review_stat_gds(n = 10L, na_first = FALSE, row_data = rd)
  pv <- preview(g, n = 5)
  expect_equal(pv$space$labels, paste0("r", 1:5))
  expect_equal(nrow(row_data(pv)), 5L)
  expect_equal(row_data(pv)$roi, paste0("roi", 1:5))
  expect_equal(dim(assay(pv, "beta"))[1L], 5L)
})

test_that("new_gds validates sample label length", {
  expect_error(
    new_gds(list(beta = array(1, c(2, 1, 1)), var = array(1, c(2, 1, 1))), space_sample_labels(c("a", "b", "c")), "s1", "c1"),
    "labels must match samples"
  )
  expect_error(
    new_gds(list(beta = array(1, c(2, 1, 1)), var = array(1, c(2, 1, 1))), space_parcels("a"), "s1", "c1"),
    "labels must match samples"
  )
})

# 5. save_plan()/load_plan() -------------------------------------------------

.review_csv <- function() {
  tmp <- tempfile(fileext = ".csv")
  rows <- c("sample,subject,contrast,beta,var")
  set.seed(7)
  for (r in c("R1", "R2", "R3")) for (s in c("s1", "s2", "s3")) for (cn in c("a", "b", "c")) {
    rows <- c(rows, sprintf("%s,%s,%s,%.4f,%.4f", r, s, cn, stats::rnorm(1), stats::runif(1, 0.1, 0.3)))
  }
  writeLines(rows, tmp)
  tmp
}

test_that("subset + reduce plans round-trip through save_plan/load_plan", {
  csv <- .review_csv()
  path <- tempfile(fileext = ".json")
  on.exit(unlink(c(csv, path)), add = TRUE)

  plan <- gds(csv) |> subset(contrast = c("a", "b"), subject = c("s1", "s3")) |> reduce("fixed")
  save_plan(plan, path)
  expect_silent(loaded <- load_plan(path))
  expect_identical(loaded$nodes[[1]]$contrast, c("a", "b"))
  expect_identical(loaded$nodes[[1]]$subject, c("s1", "s3"))
  a <- compute(plan)
  b <- compute(loaded)
  expect_equal(assays(b), assays(a))
  expect_equal(contrasts(b), c("a", "b"))
})

test_that("plans with register_map/align/map_to (incl. Matrix maps) round-trip", {
  csv <- .review_csv()
  path <- tempfile(fileext = ".json")
  on.exit(unlink(c(csv, path)), add = TRUE)

  base <- gds(csv)
  src_space <- space(base)
  ops <- setNames(rep(list(matrix(c(1, 0, 0, 0, 0.5, 0.5), nrow = 2, byrow = TRUE)), 3), c("s1", "s2", "s3"))
  fam <- MapFamily("to_two", src_space, space_sample_labels(c("T1", "T2")), by_subject = ops)
  M <- Matrix::Matrix(matrix(c(0.5, 0.5), nrow = 1), sparse = TRUE)
  plan <- base |>
    register_map(fam) |>
    align("to_two") |>
    map_to(space_sample_labels("G"), map = M) |>
    reduce("fixed")
  save_plan(plan, path)
  expect_silent(loaded <- load_plan(path))
  expect_s3_class(loaded$meta$map_families$to_two, "gds_map_family")
  expect_equal(assays(compute(loaded)), assays(compute(plan)))

  lin <- map_linear(space_sample_labels(c("T1", "T2")), space_sample_labels("G"), matrix(c(0.5, 0.5), nrow = 1))
  plan2 <- base |> align(fam) |> map_to(space_sample_labels("G"), map = lin)
  save_plan(plan2, path)
  loaded2 <- load_plan(path)
  expect_s3_class(loaded2$nodes[[2]]$map, "map_linear")
  expect_equal(assays(compute(loaded2)), assays(compute(plan2)))
})

test_that("save_plan errors on custom mask functions and load_plan verifies the digest", {
  csv <- .review_csv()
  path <- tempfile(fileext = ".json")
  on.exit(unlink(c(csv, path)), add = TRUE)

  custom <- gds(csv) |> mask(MaskPolicy(rule = "custom", custom = function(a) rep(TRUE, dim(a[[1]])[1])))
  expect_error(save_plan(custom, path), "custom")

  save_plan(gds(csv) |> reduce("fixed"), path)
  txt <- readLines(path)
  txt <- sub("\"method\": \"fixed\"", "\"method\": \"random\"", txt, fixed = TRUE)
  writeLines(txt, path)
  expect_warning(load_plan(path), "digest mismatch")
})

# 6. subset on dense voxel spaces --------------------------------------------

test_that("subset() on dense voxel GDS: subject subsets don't warn, sample subsets compute", {
  sp <- space_voxel(dim = c(2L, 2L, 2L), affine = diag(4))
  g <- .review_stat_gds(n = 8L, na_first = FALSE, space = sp)
  expect_no_warning(out <- compute(subset(as_plan(g), subject = "s1")))
  expect_equal(subjects(out), "s1")
  expect_identical(out$space$storage, "dense")

  out2 <- compute(subset(as_plan(g), sample = 1:4))
  expect_equal(dim(assay(out2, "beta"))[1L], 4L)
  expect_identical(out2$space$storage, "packed")
  expect_equal(out2$space$mask_idx, 1:4)
})

# 7. subset() on gds and gds_source -------------------------------------------

test_that("subset() works on realised gds and gds_source objects", {
  g <- .review_stat_gds(na_first = FALSE)
  p <- subset(g, subject = c("s1", "s2"))
  expect_s3_class(p, "gds_plan")
  expect_equal(subjects(compute(p)), c("s1", "s2"))

  src <- as_plan(g)$source
  expect_s3_class(src, "gds_source")
  p2 <- subset(src, sample = 2:3)
  expect_s3_class(p2, "gds_plan")
  expect_equal(compute(p2)$space$labels, c("r2", "r3"))
})

# 8. subjects()/contrasts() on plans ------------------------------------------

test_that("subjects()/contrasts() on plans reflect subset and reduce nodes", {
  csv <- .review_csv()
  on.exit(unlink(csv), add = TRUE)
  p <- gds(csv) |> subset(subject = c("s3", "s1"), contrast = "b")
  expect_equal(subjects(p), subjects(compute(p)))
  expect_equal(subjects(p), c("s3", "s1"))
  expect_equal(contrasts(p), contrasts(compute(p)))
  r <- reduce(p, "fixed")
  expect_equal(subjects(r), subjects(compute(r)))
  expect_equal(subjects(r), "meta")
  expect_equal(contrasts(r), "b")
})

# 9. Lancaster per-sample df ---------------------------------------------------

test_that("combine:lancaster uses per-sample degrees of freedom", {
  p <- array(c(0.01, 0.20, 0.03, 0.04, 0.50, 0.02), dim = c(2, 3, 1))
  df <- array(c(5, 40, 5, 40, 5, 40), dim = c(2, 3, 1))
  g <- new_gds(list(p = p, df = df), space_sample_labels(c("a", "b")), c("s1", "s2", "s3"), "c1")
  out <- compute(reduce(g, "combine:lancaster"))
  manual <- vapply(1:2, function(i) {
    w <- df[i, , 1]
    x2 <- sum(stats::qchisq(1 - p[i, , 1], df = 2 * w))
    stats::pchisq(x2, 2 * sum(w), lower.tail = FALSE)
  }, numeric(1))
  expect_equal(assay(out, "p_g")[, 1, 1], manual)
})

test_that("fdr:spatial grouping matches a direct Simes + weighted BH computation", {
  p <- array(c(0.001, 0.02, 0.3, 0.04, 0.5, 0.9), dim = c(6, 1, 1))
  g <- new_gds(list(p = p), space_sample_labels(paste0("v", 1:6)), "meta", "c1")
  grp <- c(1, 1, 2, 2, 3, NA)
  q <- assay(compute(posthoc(g, "fdr:spatial", options = list(group = grp))), "q")[, 1, 1]
  simes <- function(pv) min(sort(pv) * length(pv) / seq_along(pv))
  pg <- c(simes(p[1:2]), simes(p[3:4]), simes(p[5]))
  w <- c(2, 2, 1) / mean(c(2, 2, 1))
  ord <- order(pg / w)
  adj <- (3 * (pg / w)[ord]) / seq_along(ord)
  adj <- pmin(rev(cummin(rev(adj))), 1)
  qg <- numeric(3); qg[ord] <- adj
  expect_equal(q, c(qg[1], qg[1], qg[2], qg[2], qg[3], NA))
})

# 10. Dead optimizer helpers -----------------------------------------------------

test_that("dead optimizer subset helpers are removed", {
  ns <- asNamespace("fmrigds")
  expect_false(exists(".combine_subsets", envir = ns, inherits = FALSE))
  expect_false(exists(".merge_subset", envir = ns, inherits = FALSE))
})

test_that("contrasts() falls back to stats::contrasts for factors", {
  f <- factor(c("a", "b", "c"))
  expect_equal(contrasts(f), stats::contrasts(f))
  expect_equal(contrasts(f, contrasts = FALSE), stats::contrasts(f, contrasts = FALSE))
})

test_that("internal plan plumbing is not exported", {
  exports <- getNamespaceExports("fmrigds")
  internal <- c("plan", "space_voxels", "add_op", "op_reduce", "digest_plan",
                "gds_plan", "canonicalize_node", "reduce_eager", "subset_eager")
  expect_false(any(internal %in% exports))
})

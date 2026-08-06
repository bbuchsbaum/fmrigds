test_that("optimizer preserves operation order", {
  nodes <- list(
    op_map(space_parcels("grp"), matrix(1), UncertaintyRule("independent")),
    op_subset_axis(sample = 1:5),
    op_derive(c("var")),
    op_mask_policy(MaskPolicy())
  )
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt <- optimize_plan(plan)
  expect_equal(vapply(opt$nodes, `[[`, character(1), "op"),
               c("map", "subset_axis", "derive", "mask_policy"))
})

test_that("optimizer preserves distinct derive nodes", {
  nodes <- list(op_derive("var"), op_derive("t"))
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt <- optimize_plan(plan)
  expect_equal(length(opt$nodes), 2)
  expect_equal(vapply(opt$nodes, function(x) x$what, character(1)), c("var", "t"))
})

test_that("optimizer does not merge or push subset nodes", {
  nodes <- list(
    op_map(space_parcels("grp"), matrix(1), UncertaintyRule("independent")),
    op_subset_axis(sample = 1:5),
    op_subset_axis(sample = 2:4, subject = "s1")
  )
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt <- optimize_plan(plan)
  expect_equal(length(opt$nodes), 3)
  expect_equal(vapply(opt$nodes, `[[`, character(1), "op"),
               c("map", "subset_axis", "subset_axis"))
})

test_that("optimizer does not move subset across alignment", {
  family <- MapFamily(
    name = "id",
    from_space = space_parcels(c("roi1", "roi2")),
    to_space = space_parcels(c("roi1", "roi2")),
    type = "linear",
    by_subject = list(s1 = diag(2))
  )
  nodes <- list(
    op_align_to_group(family),
    op_subset_axis(sample = 1:2)
  )
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt <- optimize_plan(plan)
  expect_equal(vapply(opt$nodes, `[[`, character(1), "op"), c("align_to_group", "subset_axis"))
})

test_that("sequential positional and label subsets compose on current axes", {
  beta <- array(seq_len(5 * 3), c(5, 3, 1))
  g <- new_gds(
    assays = list(beta = beta, var = array(1, dim(beta))),
    space = space_sample_labels(letters[1:5]),
    subjects = c("s1", "s2", "s3"),
    contrasts = "c1"
  )

  result <- as_plan(g) |>
    subset(sample = 2:5, subject = c("s3", "s1")) |>
    subset(sample = 2:3, subject = "s1") |>
    compute()

  expect_equal(sample_labels(result), c("c", "d"))
  expect_equal(subjects(result), "s1")
  expect_equal(drop(assay(result, "beta")), beta[3:4, 1, 1])
})

test_that("a subset after mapping addresses the transformed sample axis", {
  g <- new_gds(
    assays = list(
      beta = array(1:3, c(3, 1, 1)),
      var = array(1, c(3, 1, 1))
    ),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = "s1",
    contrasts = "c1"
  )
  map <- rbind(first = c(1, 0, 0), second = c(0, 1, 1))

  result <- as_plan(g) |>
    map_to(space_sample_labels(c("first", "second")), map) |>
    subset(sample = 2L) |>
    compute()

  expect_equal(sample_labels(result), "second")
  expect_equal(drop(assay(result, "beta")), 5)
})

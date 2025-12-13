test_that("optimizer orders nodes by priority", {
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
               c("subset_axis", "mask_policy", "derive", "map"))
})

test_that("optimizer coalesces derive nodes", {
  nodes <- list(op_derive("var"), op_derive("t"))
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt <- optimize_plan(plan)
  expect_equal(length(opt$nodes), 1)
  expect_equal(sort(opt$nodes[[1]]$what), c("t", "var"))
})

test_that("optimizer merges subset nodes and pushes them first", {
  nodes <- list(
    op_map(space_parcels("grp"), matrix(1), UncertaintyRule("independent")),
    op_subset_axis(sample = 1:5),
    op_subset_axis(sample = 2:4, subject = "s1")
  )
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt <- optimize_plan(plan)
  expect_equal(opt$nodes[[1]]$op, "subset_axis")
  expect_equal(opt$nodes[[1]]$sample, 2:4)
  expect_equal(opt$nodes[[1]]$subject, "s1")
  expect_equal(vapply(opt$nodes, `[[`, character(1), "op"), c("subset_axis", "map"))
})

test_that("optimizer keeps align after subset", {
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
  expect_equal(vapply(opt$nodes, `[[`, character(1), "op"), c("subset_axis", "align_to_group"))
})

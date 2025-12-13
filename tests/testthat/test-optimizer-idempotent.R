test_that("optimize_plan is idempotent", {
  nodes <- list(
    op_map(space_parcels("grp"), matrix(1), UncertaintyRule("independent")),
    op_subset_axis(sample = 1:5),
    op_derive(c("var")),
    op_mask_policy(MaskPolicy())
  )
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- nodes
  opt1 <- optimize_plan(plan)
  opt2 <- optimize_plan(opt1)
  can <- function(pl) lapply(pl$nodes, canonicalize_node)
  expect_equal(can(opt1), can(opt2))
})


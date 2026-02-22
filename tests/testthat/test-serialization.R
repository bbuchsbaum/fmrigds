test_that("plan serialization roundtrip", {
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  family <- MapFamily("align", space_parcels("source"), space_parcels("target"), type = "linear", by_subject = list(s1 = matrix(1)))
  plan$nodes <- list(
    op_subset_axis(sample = 1:3),
    op_derive("var"),
    op_align_to_group(family),
    op_mask_policy(MaskPolicy()),
    op_map(space_parcels("target"), matrix(1), UncertaintyRule("independent"), combine = NULL),
    op_reduce("fixed", "1/var", by = "contrast", options = list()),
    list(op = "posthoc", method = "fdr:bh")
  )

  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)

  save_plan(plan, path)
  loaded <- load_plan(path)

  expect_equal(vapply(loaded$nodes, `[[`, character(1), "op"), c("subset_axis", "derive", "align_to_group", "mask_policy", "map", "reduce", "posthoc"))
})

test_that("plan serialization preserves posthoc options", {
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  plan$nodes <- list(
    list(
      op = "posthoc",
      method = "nt:cluster_fdr_perm",
      options = list(
        n_perm = 64L,
        q = 0.05,
        cluster_thresh = 2.7,
        tail = "two",
        two_sided_policy = "BH_all"
      )
    )
  )

  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)

  save_plan(plan, path)
  loaded <- load_plan(path)

  expect_equal(length(loaded$nodes), 1L)
  expect_equal(loaded$nodes[[1]]$op, "posthoc")
  expect_equal(loaded$nodes[[1]]$method, "nt:cluster_fdr_perm")
  expect_equal(loaded$nodes[[1]]$options$n_perm, 64L)
  expect_equal(loaded$nodes[[1]]$options$q, 0.05)
  expect_equal(loaded$nodes[[1]]$options$cluster_thresh, 2.7)
  expect_equal(loaded$nodes[[1]]$options$tail, "two")
  expect_equal(loaded$nodes[[1]]$options$two_sided_policy, "BH_all")
})

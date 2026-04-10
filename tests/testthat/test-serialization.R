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

test_that("plan serialization preserves reduce options, formula, map target space, and write nodes", {
  probe <- list(
    assays = c("beta", "var"),
    dims = gds_dims(sample = 1, subject = 1, contrast = 1),
    subjects = "s1",
    contrasts = "c1",
    space = space_sample_labels("ROI_1"),
    maps = list(),
    metadata = list(),
    columns = list()
  )
  plan <- gds_plan(gds_source("memory", list(beta = array(1, c(1, 1, 1)), var = array(1, c(1, 1, 1))), probe))
  plan$nodes <- list(
    op_map(
      target_space = space_parcels(c("G1", "G2")),
      map = matrix(c(1, 0), nrow = 2),
      uncertainty = UncertaintyRule("independent", df_rule = "none"),
      combine = "fisher"
    ),
    op_reduce(
      method = "meta:fe_reg",
      weights = "1/var",
      by = "contrast",
      options = list(eps = 1e-8),
      formula = "~ age"
    ),
    op_write(path = "out.csv", format = "csv", options = list(stats = "beta"))
  )

  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)

  save_plan(plan, path)
  loaded <- load_plan(path)

  expect_equal(loaded$nodes[[1]]$combine, "fisher")
  expect_s3_class(loaded$nodes[[1]]$target_space, "space_parcels")
  expect_equal(loaded$nodes[[1]]$target_space$labels, c("G1", "G2"))
  expect_equal(loaded$nodes[[1]]$uncertainty$df_rule, "none")
  expect_equal(loaded$nodes[[2]]$options$eps, 1e-8)
  expect_equal(loaded$nodes[[2]]$formula, "~ age")
  expect_equal(loaded$nodes[[3]]$path, "out.csv")
  expect_equal(loaded$nodes[[3]]$format, "csv")
  expect_equal(loaded$nodes[[3]]$options$stats, "beta")
})

test_that("plan serialization preserves contrast metadata and remains computable", {
  tmp <- tempfile(fileext = ".csv")
  out <- tempfile(fileext = ".json")
  on.exit(unlink(c(tmp, out)), add = TRUE)

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,baseline,0.5,0.04",
    "ROI_1,sub-01,task,0.8,0.04",
    "ROI_1,sub-02,baseline,0.4,0.05",
    "ROI_1,sub-02,task,0.9,0.05"
  ), tmp)

  plan <- gds(tmp)
  plan <- with_contrast_data(
    plan,
    data.frame(
      condition = c(0, 1),
      row.names = c("baseline", "task"),
      stringsAsFactors = FALSE
    )
  )

  save_plan(plan, out)
  loaded <- load_plan(out)

  expect_equal(contrast_data(loaded), contrast_data(plan))
  gout <- compute(loaded)
  expect_equal(contrast_data(gout), contrast_data(plan))
})

test_that("plan serialization preserves subset nodes with omitted axes", {
  tmp <- tempfile(fileext = ".csv")
  plan_file <- tempfile(fileext = ".json")
  on.exit(unlink(c(tmp, plan_file)), add = TRUE)

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,baseline,0.5,0.04",
    "ROI_1,sub-01,task,0.8,0.04",
    "ROI_1,sub-02,baseline,0.4,0.05",
    "ROI_1,sub-02,task,0.9,0.05",
    "ROI_2,sub-01,baseline,0.6,0.03",
    "ROI_2,sub-01,task,0.7,0.03",
    "ROI_2,sub-02,baseline,0.55,0.02",
    "ROI_2,sub-02,task,0.95,0.02"
  ), tmp)

  plan <- gds(tmp) |>
    subset(contrast = "task") |>
    reduce(method = "fixed")

  save_plan(plan, plan_file)
  loaded <- load_plan(plan_file)

  expect_null(loaded$nodes[[1]]$sample)
  expect_null(loaded$nodes[[1]]$subject)
  expect_equal(loaded$nodes[[1]]$contrast, "task")

  preview_out <- preview(loaded, n = 1)
  expect_equal(subjects(preview_out), "meta")
  expect_equal(contrasts(preview_out), "task")
  expect_equal(dim(assay(preview_out, "beta")), c(1, 1, 1))
})

test_that("make_linear_family builds a usable family and apply_align works", {
  beta <- array(1:6, dim = c(3, 2, 1))
  var <- array(0.1, dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var)
  from <- space_parcels(paste0("r", 1:3))
  to <- space_parcels(paste0("g", 1:3))
  ops <- list(s1 = diag(3), s2 = diag(3))
  fam <- make_linear_family("lin", from, to, ops)
  res <- fmrigds:::apply_align(fam, arrays, subjects = c("s1","s2"), space = from)
  expect_equal(dim(res$arrays$beta), c(3, 2, 1))
  expect_s3_class(res$space, "gds_space")
})

test_that("make_warp_family loads operators from paths via loader", {
  # write two small matrices to disk, read back via loader
  m1 <- diag(2); m2 <- matrix(c(0,1,1,0), nrow=2)
  p1 <- tempfile(fileext = ".csv"); p2 <- tempfile(fileext = ".csv")
  on.exit(unlink(c(p1,p2)))
  write.table(m1, file = p1, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(m2, file = p2, sep = ",", row.names = FALSE, col.names = FALSE)
  loader <- function(path) as.matrix(utils::read.csv(path, header = FALSE))

  from <- space_parcels(c("a","b"))
  to <- space_parcels(c("x","y"))
  fam <- make_warp_family("warp", from, to, c(s1 = p1, s2 = p2), loader)
  beta <- array(1:4, dim = c(2,2,1))
  var <- array(0.2, dim = c(2,2,1))
  arrays <- list(beta = beta, var = var)
  res <- fmrigds:::apply_align(fam, arrays, subjects = c("s1","s2"), space = from)
  expect_equal(res$arrays$beta[,1,1], c(1,2))
  expect_equal(res$arrays$beta[,2,1], c(4,3))
})

test_that("register_alignment sugar works on a plan", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,s1,c1,1,0.1",
    "ROI_2,s1,c1,2,0.1",
    "ROI_1,s2,c1,3,0.1",
    "ROI_2,s2,c1,4,0.1"
  ), tmp)
  plan <- gds(tmp)
  from <- space_parcels(c("ROI_1","ROI_2"))
  to <- from
  ops <- list(s1 = diag(2), s2 = diag(2))
  fam <- make_linear_family("id", from, to, ops)
  plan2 <- register_alignment(plan, fam)
  expect_true("id" %in% list_alignments(plan2))
  fam2 <- get_alignment(plan2, "id")
  expect_s3_class(fam2, "gds_map_family")
})


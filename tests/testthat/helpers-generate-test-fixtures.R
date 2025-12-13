#' Generate Small Synthetic fmristore Test Fixtures
#'
#' Creates minimal HDF5 files for testing fmristore adapter integration.
#' Files are kept small (<100KB) for fast test execution.
#'
#' Run this script manually when test fixtures need to be regenerated:
#' source("tests/testthat/helpers-generate-test-fixtures.R")
#' generate_all_fmristore_fixtures()

.fixture_path <- function(path) {
  # If path already exists as-is, use it directly (e.g., when called from package root)
  if (file.exists(path) || file.exists(dirname(path))) {
    return(path)
  }
  # Otherwise use test_path for running within testthat context
  if (requireNamespace("testthat", quietly = TRUE)) {
    testthat::test_path(path)
  } else {
    path
  }
}

#' Generate LabeledVolumeSet test fixture
#' 
#' @param output_path Path to output HDF5 file
#' @param dims Spatial dimensions (default: c(4, 5, 3))
#' @param n_contrasts Number of contrasts (default: 2)
#' @param labels Contrast labels (default: c("contrast_A", "contrast_B"))
#' @return Path to created file (invisible)
generate_labeled_volume_fixture <- function(
  output_path = "tests/testthat/testdata/labeled_volume_small.h5",
  dims = c(4, 5, 3),
  n_contrasts = 2,
  labels = c("contrast_A", "contrast_B")
) {
  output_path <- .fixture_path(output_path)
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package required")
  }
  if (!requireNamespace("fmristore", quietly = TRUE)) {
    stop("fmristore package required")
  }

  message("Generating LabeledVolumeSet fixture: ", output_path)

  # Create 4D data array [X, Y, Z, nVols]
  set.seed(12345)  # Reproducible test data
  arr_data <- array(rnorm(prod(dims) * n_contrasts),
                   dim = c(dims, n_contrasts))

  # Create NeuroSpace
  spc <- neuroim2::NeuroSpace(
    dim = c(dims, n_contrasts),
    spacing = c(2.0, 2.0, 2.0)  # 2mm isotropic
  )

  # Create DenseNeuroVec
  vec_obj <- neuroim2::DenseNeuroVec(arr_data, space = spc)

  # Create mask (60% valid voxels)
  mask_array <- array(
    sample(c(TRUE, FALSE), prod(dims),
           replace = TRUE, prob = c(0.6, 0.4)),
    dim = dims
  )
  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_array,
    space = neuroim2::drop_dim(spc)
  )

  # Write to HDF5
  h5_handle <- fmristore::write_labeled_vec(
    vec = vec_obj,
    mask = mask_vol,
    labels = labels,
    file = output_path,
    compression = 4
  )

  h5_handle$close_all()

  file_size <- file.info(output_path)$size
  message("  Created: ", basename(output_path), " (",
          round(file_size / 1024, 1), " KB)")

  return(invisible(output_path))
}


#' Generate H5ParcellatedScan test fixture
#'
#' @param output_path Path to output HDF5 file
#' @param mask_dims Spatial dimensions (default: c(4, 5, 3))
#' @param n_clusters Number of parcels/ROIs (default: 3)
#' @param n_time_run1 Timepoints for run 1 (default: 10)
#' @param n_time_run2 Timepoints for run 2 (default: 12)
#' @return Path to created file (invisible)
generate_parcellated_fixture <- function(
  output_path = "tests/testthat/testdata/parcellated_small.h5",
  mask_dims = c(4, 5, 3),
  n_clusters = 3,
  n_time_run1 = 10,
  n_time_run2 = 12
) {
  output_path <- .fixture_path(output_path)
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package required")
  }
  if (!requireNamespace("fmristore", quietly = TRUE)) {
    stop("fmristore package required")
  }

  message("Generating H5ParcellatedScan fixture: ", output_path)

  # Use fmristore's existing comprehensive helper
  # This creates a complete multi-scan parcellated file
  fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(
    file_path = output_path,
    master_mask_dims = mask_dims,
    num_master_clusters = n_clusters,
    n_time_run1 = n_time_run1,
    n_time_run2 = n_time_run2
  )

  # WORKAROUND: fmristore helper writes cluster_map as 1D array of IDs
  # We need to convert it to 3D cluster_map matching mask dimensions
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package required")
  }

  h5 <- hdf5r::H5File$new(output_path, mode = "r+")
  on.exit(h5$close_all())

  # Read existing 1D cluster_ids
  if (h5$exists("/cluster_ids")) {
    cluster_ids <- as.integer(h5$open("/cluster_ids")$read())

    # Read mask to get dimensions
    mask <- h5$open("/mask")$read()
    mask_dims_actual <- dim(mask)

    # Create 3D cluster_map by assigning voxels to clusters
    set.seed(12345)  # Reproducible
    cluster_map_3d <- array(0L, dim = mask_dims_actual)

    # Assign masked voxels to clusters
    masked_idx <- which(mask > 0)
    n_masked <- length(masked_idx)

    if (n_masked > 0 && length(cluster_ids) > 0) {
      # Distribute masked voxels among clusters
      assignments <- sample(cluster_ids, n_masked, replace = TRUE)
      cluster_map_3d[masked_idx] <- assignments
    }

    # Delete old 1D cluster_map if exists
    if (h5$exists("/cluster_map")) {
      h5$link_delete("/cluster_map")
    }

    # Write new 3D cluster_map
    # Ensure cluster_map_3d has correct dimensions and storage mode
    cluster_map_3d <- array(as.integer(cluster_map_3d), dim = mask_dims_actual)

    # Use simple bracket assignment which should preserve dimensions
    h5[["cluster_map"]] <- cluster_map_3d

    message("  Fixed: cluster_map is now 3D ", paste(mask_dims_actual, collapse="x"))
  }

  # Explicitly flush and close
  if (h5$is_valid) {
    h5$flush()
    h5$close_all()
  }
  on.exit(NULL)

  file_size <- file.info(output_path)$size
  message("  Created: ", basename(output_path), " (",
          round(file_size / 1024, 1), " KB)")

  return(invisible(output_path))
}


#' Generate LatentNeuroVec test fixture
#'
#' @param output_path Path to output HDF5 file
#' @param dims Spatial + temporal dimensions (default: c(4, 5, 3, 10))
#' @param k Number of basis components (default: 3)
#' @return Path to created file (invisible)
generate_latent_vec_fixture <- function(
  output_path = "tests/testthat/testdata/latent_vec_small.lv.h5",
  dims = c(4, 5, 3, 10),  # unused for minimal file
  k = 3                    # number of basis components
) {
  output_path <- .fixture_path(output_path)
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package required")
  }
  set.seed(12345)
  V <- 5L
  basis_matrix <- matrix(rnorm(k * V), nrow = k, ncol = V)
  h5 <- hdf5r::H5File$new(output_path, mode = "w")
  on.exit({ if (h5$is_valid) h5$close_all() }, add = TRUE)
  basis_grp <- h5$create_group("basis")
  basis_grp$create_dataset("basis_matrix", basis_matrix)
  basis_grp$create_dataset("basis_method", "ICA")
  scans_grp <- h5$create_group("scans")
  scans_grp$create_group("sub-01")$create_dataset("embedding", matrix(rnorm(k), nrow = k, ncol = 1))
  scans_grp$create_group("sub-02")$create_dataset("embedding", matrix(rnorm(k), nrow = k, ncol = 1))
  scans_grp$close(); basis_grp$close(); h5$close_all()
  invisible(output_path)
}


#' Generate all fmristore test fixtures
#'
#' Creates all three layout types in tests/testthat/testdata/
#'
#' @return List of created file paths (invisible)
#' @export
generate_all_fmristore_fixtures <- function() {
  message("\n=== Generating fmristore Test Fixtures ===\n")

  # Ensure testdata directory exists
  testdata_dir <- "tests/testthat/testdata"
  if (!dir.exists(.fixture_path(testdata_dir))) {
    dir.create(.fixture_path(testdata_dir), recursive = TRUE)
    message("Created directory: ", testdata_dir, "\n")
  }

  fixtures <- list()

  # 1. LabeledVolumeSet
  tryCatch({
    fixtures$labeled_volume <- generate_labeled_volume_fixture()
  }, error = function(e) {
    warning("Failed to generate LabeledVolumeSet fixture: ", e$message)
  })

  # 2. H5ParcellatedScan
  tryCatch({
    fixtures$parcellated <- generate_parcellated_fixture()
  }, error = function(e) {
    warning("Failed to generate H5ParcellatedScan fixture: ", e$message)
  })

  # 3. LatentNeuroVec
  tryCatch({
    fixtures$latent_vec <- generate_latent_vec_fixture()
  }, error = function(e) {
    warning("Failed to generate LatentNeuroVec fixture: ", e$message)
  })

  message("\n=== Fixture Generation Complete ===")
  message("Total fixtures created: ", length(fixtures))

  total_size <- sum(vapply(fixtures, function(f) {
    if (file.exists(f)) file.info(f)$size else 0
  }, numeric(1)))
  message("Total size: ", round(total_size / 1024, 1), " KB\n")

  return(invisible(fixtures))
}


# If run directly, generate all fixtures
if (interactive()) {
  message("To generate test fixtures, run:")
  message("  generate_all_fmristore_fixtures()")
}

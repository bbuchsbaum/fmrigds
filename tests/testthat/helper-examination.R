.group_examination_fixture <- function() {
  n_sample <- 40L
  n_subject <- 10L
  signal <- sin(seq(0, 2 * pi, length.out = n_sample))
  beta <- array(NA_real_, c(n_sample, n_subject, 1L))
  for (i in 1:8) {
    beta[, i, 1] <- signal + 0.03 * cos(seq_len(n_sample) + i)
  }
  beta[, 9, 1] <- signal + 3 * sin(seq_len(n_sample) * 1.7)
  beta[, 10, 1] <- -signal
  var <- array(0.04, dim(beta))
  var[, 8, 1] <- 0.001
  var[, 9, 1] <- 9
  new_gds(
    list(beta = beta, var = var),
    space_sample_labels(paste0("feature-", seq_len(n_sample))),
    paste0("s", seq_len(n_subject)),
    "task"
  )
}

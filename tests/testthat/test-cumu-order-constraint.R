## The reparameterisation alpha_t = atilde_1 - sum_{s=2..t} exp(atilde_s)
## guarantees alpha_1 >= alpha_2 >= ... >= alpha_T by construction. This
## test guards against a back-transform bug.

test_that("fitted alpha respects the proportional-odds ordering", {
  skip_on_cran()
  n_rep <- 15L
  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = 250L, p = 2L, T = 2L,
                            alpha = c(0.8, -0.4), beta = c(0.4, -0.3),
                            seed = 500L + r)
    M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                   q = fx$q, capture = "B", seed = 600L + r)
    theta_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
    res <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                 theta_estimated = theta_init,
                                 q_vec = fx$q, T = 2L)
    if (res[length(res)] != 1L) next
    alpha_hat <- PACS:::alpha_from_atilde(res[1:2])
    expect_gte(alpha_hat[1], alpha_hat[2])
  }
})


test_that("alpha_from_atilde / atilde_jacobian are consistent with FD", {
  set.seed(7L)
  for (T in 1:4) {
    atilde <- rnorm(T, sd = 0.3)
    J <- PACS:::atilde_jacobian(atilde)
    h <- 1e-6
    J_fd <- matrix(0, nrow = T, ncol = T)
    for (s in seq_len(T)) {
      a_p <- atilde; a_p[s] <- a_p[s] + h
      a_m <- atilde; a_m[s] <- a_m[s] - h
      J_fd[, s] <- (PACS:::alpha_from_atilde(a_p) -
                    PACS:::alpha_from_atilde(a_m)) / (2 * h)
    }
    expect_equal(J, J_fd, tolerance = 1e-6)
  }
})

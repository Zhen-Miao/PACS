# Cross-validation: generate test data, run R PACS, save results
# This script sources the R functions directly to avoid needing install

# Load required packages
library(parallel)

# Source all R functions
r_dir <- "/home/user/PACS/R"
source(file.path(r_dir, "param-estimate_logit_get_p_by_t_June.R"))
source(file.path(r_dir, "differential_identification.R"))
source(file.path(r_dir, "PACS_test_logit.R"))
source(file.path(r_dir, "PACS_test_cumulative.R"))
source(file.path(r_dir, "PACS_test_auto.R"))
source(file.path(r_dir, "cell_type_annotation.R"))
source(file.path(r_dir, "CCT_function.R"))

library(Matrix)

out_dir <- "/home/user/PACS/pacs-py/tests/cross_validation"

# ======================================================================
# 1. Generate simple test data (reproducible)
# ======================================================================
set.seed(42)
n_cells <- 60
n_features <- 10

cell_types <- c(rep("A", 30), rep("B", 30))
cap_rates <- runif(n_cells, 0.3, 0.7)

# PIC matrix: first 5 features are DARs
pic_matrix <- matrix(0, nrow = n_features, ncol = n_cells)
rownames(pic_matrix) <- paste0("f_", 1:n_features)

for (i in 1:5) {
  p_A <- 0.6 + i * 0.05  # 0.65, 0.70, 0.75, 0.80, 0.85
  p_B <- 0.05 + i * 0.02 # 0.07, 0.09, 0.11, 0.13, 0.15
  for (j in 1:n_cells) {
    p <- ifelse(cell_types[j] == "A", p_A, p_B)
    pic_matrix[i, j] <- ifelse(runif(1) < p * cap_rates[j], 1, 0)
  }
}
for (i in 6:n_features) {
  p_both <- 0.2 + i * 0.02
  for (j in 1:n_cells) {
    pic_matrix[i, j] <- ifelse(runif(1) < p_both * cap_rates[j], 1, 0)
  }
}

# Save test data
write.csv(pic_matrix, file.path(out_dir, "pic_matrix.csv"))
write.csv(data.frame(cell_type = cell_types), file.path(out_dir, "metadata.csv"),
          row.names = FALSE)
write.csv(data.frame(cap_rates = cap_rates), file.path(out_dir, "cap_rates.csv"),
          row.names = FALSE)

cat("Test data generated:\n")
cat("  pic_matrix:", dim(pic_matrix), "\n")
cat("  n_cells:", n_cells, "\n")
cat("  n_features:", n_features, "\n")

# ======================================================================
# 2. Test pacs_test_logit
# ======================================================================
cat("\n--- Running pacs_test_logit ---\n")

metadata <- data.frame(cell_type = factor(cell_types))

result_logit <- pacs_test_logit(
  covariate_meta.data = metadata,
  formula_full = ~ cell_type,
  formula_null = ~ 1,
  pic_matrix = pic_matrix,
  cap_rates = cap_rates,
  n_cores = 1
)

cat("P-values (logit):\n")
print(round(result_logit$pacs_p_val, 10))
cat("Convergence:\n")
print(result_logit$pacs_converged)

write.csv(data.frame(
  feature = names(result_logit$pacs_p_val),
  p_value = result_logit$pacs_p_val
), file.path(out_dir, "r_logit_pvalues.csv"), row.names = FALSE)

conv_mat <- matrix(result_logit$pacs_converged, ncol = 2)
write.csv(data.frame(
  null_conv = conv_mat[, 1],
  full_conv = conv_mat[, 2]
), file.path(out_dir, "r_logit_convergence.csv"), row.names = FALSE)

# ======================================================================
# 3. Test individual functions (loss, gradient, info matrix)
# ======================================================================
cat("\n--- Testing individual functions ---\n")

X <- model.matrix(~ cell_type, data = metadata)
theta <- c(0.1, 0.3)
x_times_theta <- X %*% theta
p_bg <- as.vector(1 - 1 / (exp(x_times_theta) + 1))
q_vec <- cap_rates
y_vec <- pic_matrix[1, ]

# loss_fun
lf <- loss_fun(p_bg, q_vec, y_vec)
cat("loss_fun:", lf, "\n")

# loss_gradient
lg <- loss_gradient(X, p_bg, q_vec, y_vec)
cat("loss_gradient:", lg, "\n")

# infor_mat
im <- infor_mat(X, p_bg, q_vec)
cat("infor_mat:\n")
print(im)

# loss_fun_star
lfs <- loss_fun_star(X, p_bg, q_vec, y_vec)
cat("loss_fun_star:", lfs, "\n")

# Save individual function results
results <- list(
  loss_fun = lf,
  loss_gradient = as.vector(lg),
  infor_mat = as.vector(im),
  loss_fun_star = lfs,
  theta = theta,
  p_bg = p_bg,
  y_vec = y_vec
)
saveRDS(results, file.path(out_dir, "r_individual_results.rds"))

# Also save as CSV for Python to read
write.csv(data.frame(
  loss_fun = lf,
  loss_fun_star = lfs,
  loss_gradient_1 = lg[1],
  loss_gradient_2 = lg[2],
  infor_mat_11 = im[1,1],
  infor_mat_12 = im[1,2],
  infor_mat_21 = im[2,1],
  infor_mat_22 = im[2,2]
), file.path(out_dir, "r_function_results.csv"), row.names = FALSE)

# Save the design matrix for Python to use
write.csv(X, file.path(out_dir, "design_matrix.csv"), row.names = FALSE)

# ======================================================================
# 4. Test CCT
# ======================================================================
cat("\n--- Testing CCT ---\n")

pval_mat <- matrix(c(0.01, 0.05, 0.1,
                      0.5,  0.6,  0.7,
                      1e-8, 0.3,  0.9), nrow = 3, byrow = TRUE)
rownames(pval_mat) <- paste0("f_", 1:3)
colnames(pval_mat) <- paste0("p_", 1:3)

cct_result <- CCT_internal_horizontal(pval_mat)
cat("CCT p-values:\n")
print(cct_result)

write.csv(data.frame(
  feature = names(cct_result),
  combined_p = cct_result
), file.path(out_dir, "r_cct_results.csv"), row.names = FALSE)

# Save CCT input
write.csv(pval_mat, file.path(out_dir, "cct_input.csv"))

# ======================================================================
# 5. Test cell type annotation
# ======================================================================
cat("\n--- Testing cell type annotation ---\n")

set.seed(123)
n_regions <- 30
n_types <- 2

r_by_t <- matrix(0, nrow = n_regions, ncol = n_types)
colnames(r_by_t) <- c("TypeA", "TypeB")
r_by_t[1:15, 1] <- runif(15, 0.4, 0.8)
r_by_t[16:30, 2] <- runif(15, 0.4, 0.8)
r_by_t[1:15, 2] <- runif(15, 0.02, 0.1)
r_by_t[16:30, 1] <- runif(15, 0.02, 0.1)

# Create 6 test cells (3 type A-like, 3 type B-like)
set.seed(456)
test_cells <- matrix(0, nrow = n_regions, ncol = 6)
for (j in 1:3) {
  test_cells[, j] <- rbinom(n_regions, 1, ifelse(1:n_regions <= 15, 0.7, 0.1))
}
for (j in 4:6) {
  test_cells[, j] <- rbinom(n_regions, 1, ifelse(1:n_regions <= 15, 0.1, 0.7))
}

# Run annotation manually (the R code has a bug in do.call(rbind(...)))
# We replicate the core logic here for cross-validation

r_by_t_mat <- as.matrix(r_by_t)
n_cells_annot <- ncol(test_cells)
n_reads <- colSums(test_cells)
sum_prob <- colSums(r_by_t_mat)
capturing_rate_mat <- t(outer(n_reads, sum_prob, "/"))
capturing_rate_mat[capturing_rate_mat > 0.9995] <- 0.9995
capturing_rate_mat[capturing_rate_mat < 0.00005] <- 0.00005

annot_result <- matrix(NA, nrow = n_cells_annot, ncol = ncol(r_by_t_mat))
colnames(annot_result) <- colnames(r_by_t_mat)

test_cells_mat <- as.matrix(test_cells)
for (i_cell in 1:n_cells_annot) {
  pqbyt <- r_by_t_mat * capturing_rate_mat[, i_cell]
  lg_pqbyt <- log(pqbyt)
  lg_pqbyt_q <- log1p(-1 * pqbyt)
  xvec <- test_cells_mat[, i_cell]
  annot_result[i_cell, ] <- colSums(lg_pqbyt * xvec + lg_pqbyt_q * (1 - xvec) * 1)
}

cat("Annotation log-likelihoods:\n")
print(round(annot_result, 4))

write.csv(r_by_t, file.path(out_dir, "r_by_t.csv"))
write.csv(test_cells, file.path(out_dir, "test_cells.csv"))
write.csv(annot_result, file.path(out_dir, "r_annotation_results.csv"))

cat("\n=== R cross-validation complete ===\n")

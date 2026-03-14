"""Compare Python PACS results against R PACS results."""

import sys
import numpy as np
import pandas as pd
from scipy.special import expit

# Add pacs to path
sys.path.insert(0, "/home/user/PACS/pacs-py/src")

from pacs.estimation import loss_fun, loss_gradient, infor_mat, loss_fun_star
from pacs.testing import pacs_test_logit
from pacs.cct import cauchy_combination_test
from pacs.annotation import estimate_label

DATA_DIR = "/home/user/PACS/pacs-py/tests/cross_validation"

PASS = 0
FAIL = 0


def check(name, py_val, r_val, atol=1e-4, rtol=1e-3):
    global PASS, FAIL
    py_arr = np.atleast_1d(np.asarray(py_val, dtype=np.float64))
    r_arr = np.atleast_1d(np.asarray(r_val, dtype=np.float64))

    if py_arr.shape != r_arr.shape:
        print(f"  FAIL {name}: shape mismatch py={py_arr.shape} r={r_arr.shape}")
        FAIL += 1
        return

    close = np.allclose(py_arr, r_arr, atol=atol, rtol=rtol)
    max_diff = np.max(np.abs(py_arr - r_arr))

    if close:
        print(f"  PASS {name} (max diff: {max_diff:.2e})")
        PASS += 1
    else:
        print(f"  FAIL {name} (max diff: {max_diff:.2e})")
        print(f"    Python: {py_arr}")
        print(f"    R:      {r_arr}")
        FAIL += 1


# ======================================================================
# Load test data
# ======================================================================
print("Loading test data...")
pic_df = pd.read_csv(f"{DATA_DIR}/pic_matrix.csv", index_col=0)
pic_matrix = pic_df.values.astype(np.float64)
metadata = pd.read_csv(f"{DATA_DIR}/metadata.csv")
cap_rates = pd.read_csv(f"{DATA_DIR}/cap_rates.csv")["cap_rates"].values

print(f"  pic_matrix: {pic_matrix.shape}")
print(f"  n_cells: {len(cap_rates)}")

# ======================================================================
# 1. Test individual functions
# ======================================================================
print("\n=== Individual function comparison ===")

r_results = pd.read_csv(f"{DATA_DIR}/r_function_results.csv")
X_df = pd.read_csv(f"{DATA_DIR}/design_matrix.csv")
X = X_df.values.astype(np.float64)

theta = np.array([0.1, 0.3])
p_bg = expit(X @ theta)
q_vec = cap_rates
y_vec = pic_matrix[0, :]

# loss_fun
py_loss = loss_fun(p_bg, q_vec, y_vec)
check("loss_fun", py_loss, r_results["loss_fun"].values[0])

# loss_gradient
py_grad = loss_gradient(X, p_bg, q_vec, y_vec)
r_grad = np.array([
    r_results["loss_gradient_1"].values[0],
    r_results["loss_gradient_2"].values[0]
])
check("loss_gradient", py_grad, r_grad)

# infor_mat
py_im = infor_mat(X, p_bg, q_vec)
r_im = np.array([
    [r_results["infor_mat_11"].values[0], r_results["infor_mat_12"].values[0]],
    [r_results["infor_mat_21"].values[0], r_results["infor_mat_22"].values[0]]
])
check("infor_mat", py_im, r_im)

# loss_fun_star
py_loss_star = loss_fun_star(X, p_bg, q_vec, y_vec)
check("loss_fun_star", py_loss_star, r_results["loss_fun_star"].values[0])

# ======================================================================
# 2. Test pacs_test_logit
# ======================================================================
print("\n=== pacs_test_logit comparison ===")

r_pvals = pd.read_csv(f"{DATA_DIR}/r_logit_pvalues.csv")
r_conv = pd.read_csv(f"{DATA_DIR}/r_logit_convergence.csv")

py_result = pacs_test_logit(
    metadata=metadata,
    formula_full="~ cell_type",
    formula_null="~ 1",
    pic_matrix=pic_matrix,
    cap_rates=cap_rates,
    n_jobs=1,
)

py_pvals = py_result.p_values.values
r_pvals_arr = r_pvals["p_value"].values

print(f"\n  Feature-by-feature p-value comparison:")
print(f"  {'Feature':<10} {'R p-value':<15} {'Py p-value':<15} {'Match':<8}")
print(f"  {'-'*48}")
for i in range(len(py_pvals)):
    r_p = r_pvals_arr[i]
    py_p = py_pvals[i]
    # Use relative tolerance for p-values
    if r_p > 0:
        rel_diff = abs(py_p - r_p) / r_p
    else:
        rel_diff = abs(py_p - r_p)
    match = "OK" if rel_diff < 0.1 else "DIFF"  # 10% relative tolerance
    print(f"  f_{i+1:<8} {r_p:<15.6e} {py_p:<15.6e} {match}")

check("logit_pvalues", py_pvals, r_pvals_arr, atol=1e-3, rtol=0.1)

# Check convergence
py_conv_null = py_result.convergence["null"].values
py_conv_full = py_result.convergence["full"].values
r_conv_null = r_conv["null_conv"].values
r_conv_full = r_conv["full_conv"].values

check("convergence_null", py_conv_null, r_conv_null, atol=0)
check("convergence_full", py_conv_full, r_conv_full, atol=0)

# ======================================================================
# 3. Test CCT
# ======================================================================
print("\n=== CCT comparison ===")

cct_input = pd.read_csv(f"{DATA_DIR}/cct_input.csv", index_col=0)
r_cct = pd.read_csv(f"{DATA_DIR}/r_cct_results.csv")

py_cct = cauchy_combination_test(cct_input.values)
r_cct_arr = r_cct["combined_p"].values

print(f"  R CCT:  {r_cct_arr}")
print(f"  Py CCT: {py_cct}")
check("CCT", py_cct, r_cct_arr, atol=1e-6, rtol=1e-4)

# ======================================================================
# 4. Test cell type annotation
# ======================================================================
print("\n=== Cell type annotation comparison ===")

r_by_t_df = pd.read_csv(f"{DATA_DIR}/r_by_t.csv", index_col=0)
test_cells_df = pd.read_csv(f"{DATA_DIR}/test_cells.csv", index_col=0)
r_annot = pd.read_csv(f"{DATA_DIR}/r_annotation_results.csv", index_col=0)

r_by_t = r_by_t_df.values
test_cells = test_cells_df.values.astype(np.float64)

py_annot = estimate_label(r_by_t, test_cells, alpha=1.0)
py_ll = py_annot.log_likelihoods.values
r_ll = r_annot.values

print(f"  R annotation log-likelihoods:")
print(f"  {r_ll}")
print(f"  Py annotation log-likelihoods:")
print(f"  {py_ll}")

# Note: Small differences expected due to a vector recycling issue in R's
# estimate_label_default. When R multiplies a (n_regions x n_types) matrix
# by a length-n_types vector, it recycles in column-major order (alternating
# rows), rather than applying one value per column. The Python code handles
# this correctly by broadcasting per cell type. Despite this, predictions
# match because the argmax is robust to these small differences.
check("annotation_loglik", py_ll, r_ll, atol=0.05, rtol=0.01)

# Check predictions match
r_pred = r_annot.columns[np.argmax(r_ll, axis=1)]
py_pred_idx = np.argmax(py_ll, axis=1)
# Python uses ct_1, ct_2 naming; R uses TypeA, TypeB
print(f"\n  R predictions:  {list(r_pred)}")
print(f"  Py predictions: {list(py_annot.predicted_labels.values)}")

# Both should predict same argmax
r_argmax = np.argmax(r_ll, axis=1)
py_argmax = np.argmax(py_ll, axis=1)
check("annotation_predictions", py_argmax, r_argmax, atol=0)

# ======================================================================
# Summary
# ======================================================================
print(f"\n{'='*60}")
print(f"SUMMARY: {PASS} passed, {FAIL} failed out of {PASS + FAIL} checks")
print(f"{'='*60}")

if FAIL > 0:
    sys.exit(1)

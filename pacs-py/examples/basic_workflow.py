"""Basic PACS workflow example with synthetic snATAC-seq data."""

import numpy as np
import pandas as pd

from pacs import (
    cauchy_combination_test,
    estimate_label,
    pacs_test,
)


def generate_synthetic_data(
    n_cells=100, n_features=50, n_dars=20, seed=42
):
    """Generate synthetic snATAC-seq data with known DARs."""
    rng = np.random.RandomState(seed)

    # Two cell types with equal proportions
    cell_types = np.array(
        ["TypeA"] * (n_cells // 2) + ["TypeB"] * (n_cells // 2)
    )
    metadata = pd.DataFrame({"cell_type": cell_types})

    # Capturing rates (depth effect)
    cap_rates = rng.uniform(0.3, 0.8, size=n_cells)

    # PIC matrix
    pic_matrix = np.zeros((n_features, n_cells), dtype=np.float64)

    # DARs: different open probability between cell types
    for i in range(n_dars):
        p_A = rng.uniform(0.5, 0.9)
        p_B = rng.uniform(0.01, 0.15)
        for j in range(n_cells):
            p = p_A if cell_types[j] == "TypeA" else p_B
            if rng.random() < p * cap_rates[j]:
                pic_matrix[i, j] = 1.0

    # Non-DARs: same open probability
    for i in range(n_dars, n_features):
        p = rng.uniform(0.15, 0.35)
        for j in range(n_cells):
            if rng.random() < p * cap_rates[j]:
                pic_matrix[i, j] = 1.0

    return pic_matrix, metadata, cap_rates


def main():
    print("=" * 60)
    print("PACS Python - Basic Workflow Example")
    print("=" * 60)

    # --- 1. Differential Accessible Region Testing ---
    print("\n1. Generating synthetic data...")
    pic_matrix, metadata, cap_rates = generate_synthetic_data()
    print(
        f"   Data: {pic_matrix.shape[0]} features x {pic_matrix.shape[1]} cells"
    )

    print("\n2. Running PACS test...")
    result = pacs_test(
        pic_matrix=pic_matrix,
        metadata=metadata,
        formula_full="~ cell_type",
        formula_null="~ 1",
        cap_rates=cap_rates,
        verbose=True,
    )

    print("\n3. Results summary:")
    n_significant = (result.p_values < 0.05).sum()
    print(f"   Significant features (p < 0.05): {n_significant}")
    print(f"   Top 5 p-values:")
    top5 = result.p_values.nsmallest(5)
    for name, pval in top5.items():
        print(f"     {name}: {pval:.2e}")

    # --- 2. Cell Type Annotation ---
    print("\n4. Cell type annotation example...")
    n_regions = 50
    n_types = 2

    rng = np.random.RandomState(123)
    r_by_t = np.zeros((n_regions, n_types))
    r_by_t[:25, 0] = rng.uniform(0.5, 0.8, 25)
    r_by_t[25:, 1] = rng.uniform(0.5, 0.8, 25)
    r_by_t[:25, 1] = rng.uniform(0.01, 0.1, 25)
    r_by_t[25:, 0] = rng.uniform(0.01, 0.1, 25)

    # Generate test cells
    test_cells = rng.binomial(1, 0.3, size=(n_regions, 10)).astype(np.float64)
    annotation = estimate_label(r_by_t, test_cells)
    print(f"   Predicted labels: {annotation.predicted_labels.values}")

    # --- 3. Cauchy Combination Test ---
    print("\n5. Cauchy Combination Test example...")
    p_matrix = np.array(
        [
            [0.001, 0.01, 0.1],
            [0.5, 0.6, 0.7],
            [1e-10, 0.5, 0.9],
        ]
    )
    combined = cauchy_combination_test(p_matrix)
    print(f"   Input p-values:\n{p_matrix}")
    print(f"   Combined p-values: {combined}")

    print("\n" + "=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()

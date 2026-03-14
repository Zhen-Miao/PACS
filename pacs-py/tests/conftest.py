"""Shared test fixtures for PACS tests."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def rng():
    """Seeded random number generator for reproducibility."""
    return np.random.RandomState(42)


@pytest.fixture
def small_data(rng):
    """Small synthetic snATAC-seq dataset for testing.

    Creates a dataset with 2 cell types, 50 cells, and 20 features.
    Features 0-9 are differentially accessible (DARs) and features
    10-19 are not.
    """
    n_cells = 50
    n_features = 20
    n_dars = 10

    # Cell types: 25 cells each
    cell_types = np.array(["A"] * 25 + ["B"] * 25)

    # Capturing rates
    cap_rates = rng.uniform(0.3, 0.8, size=n_cells)

    # Generate PIC matrix
    pic_matrix = np.zeros((n_features, n_cells), dtype=np.float64)

    # DARs: different open probability between cell types
    for i in range(n_dars):
        p_A = rng.uniform(0.5, 0.8)
        p_B = rng.uniform(0.05, 0.2)
        for j in range(n_cells):
            p = p_A if cell_types[j] == "A" else p_B
            pic_matrix[i, j] = 1.0 if rng.random() < p * cap_rates[j] else 0.0

    # Non-DARs: same open probability across cell types
    for i in range(n_dars, n_features):
        p_both = rng.uniform(0.2, 0.4)
        for j in range(n_cells):
            pic_matrix[i, j] = (
                1.0 if rng.random() < p_both * cap_rates[j] else 0.0
            )

    metadata = pd.DataFrame({"cell_type": cell_types})

    return {
        "pic_matrix": pic_matrix,
        "metadata": metadata,
        "cap_rates": cap_rates,
        "n_dars": n_dars,
        "n_features": n_features,
        "n_cells": n_cells,
    }


@pytest.fixture
def tiny_data():
    """Minimal 3-cell, 2-feature dataset for unit tests."""
    pic_matrix = np.array([[1, 0, 1], [0, 1, 0]], dtype=np.float64)
    cap_rates = np.array([0.5, 0.6, 0.7])
    X = np.array([[1, 0], [1, 1], [1, 0]], dtype=np.float64)
    theta = np.array([0.1, 0.2])

    return {
        "pic_matrix": pic_matrix,
        "cap_rates": cap_rates,
        "X": X,
        "theta": theta,
    }

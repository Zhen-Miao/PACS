"""Tests for the annotation module."""

import numpy as np

from pacs.annotation import estimate_label, estimate_label_selected_peaks


class TestEstimateLabel:
    def test_basic(self, rng):
        n_regions = 50
        n_cells = 20
        n_types = 3

        # Reference probability matrix
        r_by_t = rng.uniform(0.01, 0.5, size=(n_regions, n_types))

        # Binary accessibility matrix
        in_r_by_c = rng.binomial(1, 0.3, size=(n_regions, n_cells)).astype(
            np.float64
        )

        result = estimate_label(r_by_t, in_r_by_c)

        assert result.log_likelihoods.shape == (n_cells, n_types)
        assert result.predicted_labels.shape == (n_cells,)
        assert all(
            label in [f"ct_{i + 1}" for i in range(n_types)]
            for label in result.predicted_labels
        )

    def test_correct_prediction(self, rng):
        """Cells should be classified to the cell type with matching profile."""
        n_regions = 100
        n_types = 2

        # Type A has high prob in first half, type B in second half
        r_by_t = np.zeros((n_regions, n_types))
        r_by_t[:50, 0] = 0.8  # type A peaks
        r_by_t[50:, 1] = 0.8  # type B peaks
        r_by_t[:50, 1] = 0.05
        r_by_t[50:, 0] = 0.05

        # Cell that looks like type A
        cell_A = np.zeros((n_regions, 1))
        cell_A[:50, 0] = 1.0

        result = estimate_label(r_by_t, cell_A)
        assert result.predicted_labels.iloc[0] == "ct_1"

    def test_selected_peaks(self, rng):
        n_regions = 100
        n_cells = 10
        n_types = 2

        r_by_t = rng.uniform(0.01, 0.5, size=(n_regions, n_types))
        in_r_by_c = rng.binomial(1, 0.3, size=(n_regions, n_cells)).astype(
            np.float64
        )
        peaks_sel = np.array([0, 5, 10, 15, 20])

        result = estimate_label_selected_peaks(
            r_by_t, in_r_by_c, peaks_sel
        )
        assert result.log_likelihoods.shape == (n_cells, n_types)

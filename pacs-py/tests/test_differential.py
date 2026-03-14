"""Tests for the differential module."""

import numpy as np

from pacs.differential import compare_models, loss_firth_from_wii


class TestLossFirthFromWii:
    def test_basic(self):
        X = np.array([[1, 0], [1, 1], [1, 0]], dtype=np.float64)
        wii_sqrt = np.array([0.3, 0.4, 0.5])

        result = loss_firth_from_wii(wii_sqrt, X)
        assert np.isfinite(result)

    def test_positive_for_nonsingular(self):
        X = np.eye(3, dtype=np.float64)
        wii_sqrt = np.array([1.0, 1.0, 1.0])
        result = loss_firth_from_wii(wii_sqrt, X)
        # For identity matrix with unit weights, det=1, log(det)=0
        assert np.isclose(result, 0.0)


class TestCompareModels:
    def test_returns_valid_pvalues(self):
        n_cells = 30
        n_features = 3
        n_params = 2

        rng = np.random.RandomState(123)

        X_full = np.column_stack([np.ones(n_cells), rng.randn(n_cells)])
        theta_full = rng.randn(n_params, n_features) * 0.1
        theta_null = theta_full.copy()
        theta_null[1, :] = 0  # null model has second param = 0

        q_vec = rng.uniform(0.3, 0.8, n_cells)
        c_by_r = rng.binomial(1, 0.3, size=(n_cells, n_features)).astype(
            np.float64
        )

        p_values = compare_models(
            X_full=X_full,
            theta_full=theta_full,
            X_null=X_full,
            theta_null=theta_null,
            q_vec=q_vec,
            c_by_r=c_by_r,
            df_test=1,
        )

        assert p_values.shape == (n_features,)
        assert np.all((p_values >= 0) & (p_values <= 1))

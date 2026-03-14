"""Tests for the estimation module."""

import numpy as np
import pytest

from pacs.estimation import (
    ConvergenceStatus,
    _compute_wii_sqrt_X,
    infor_mat,
    irls_iter,
    irls_iter_null,
    loss_fun,
    loss_fun_star,
    loss_grad_pen,
    loss_gradient,
)


class TestLossFun:
    def test_basic(self, tiny_data):
        X = tiny_data["X"]
        theta = tiny_data["theta"]
        q_vec = tiny_data["cap_rates"]
        y_vec = tiny_data["pic_matrix"][0]

        from scipy.special import expit

        p_bg = expit(X @ theta)
        result = loss_fun(p_bg, q_vec, y_vec)
        assert np.isfinite(result)
        assert result < 0  # log-likelihood is typically negative

    def test_perfect_prediction(self):
        # When p*q matches y perfectly, loss should be highest (least negative)
        p_bg = np.array([0.99, 0.01])
        q_vec = np.array([1.0, 1.0])
        y_vec = np.array([1.0, 0.0])
        high_loss = loss_fun(p_bg, q_vec, y_vec)

        # Worse prediction
        p_bg_bad = np.array([0.5, 0.5])
        low_loss = loss_fun(p_bg_bad, q_vec, y_vec)

        assert high_loss > low_loss


class TestLossGradient:
    def test_shape(self, tiny_data):
        X = tiny_data["X"]
        theta = tiny_data["theta"]
        q_vec = tiny_data["cap_rates"]
        y_vec = tiny_data["pic_matrix"][0]

        from scipy.special import expit

        p_bg = expit(X @ theta)
        grad = loss_gradient(X, p_bg, q_vec, y_vec)
        assert grad.shape == (X.shape[1],)
        assert np.all(np.isfinite(grad))


class TestInforMat:
    def test_shape_and_symmetry(self, tiny_data):
        X = tiny_data["X"]
        theta = tiny_data["theta"]
        q_vec = tiny_data["cap_rates"]

        from scipy.special import expit

        p_bg = expit(X @ theta)
        I = infor_mat(X, p_bg, q_vec)

        assert I.shape == (X.shape[1], X.shape[1])
        np.testing.assert_allclose(I, I.T, atol=1e-12)

    def test_positive_definite(self, tiny_data):
        X = tiny_data["X"]
        theta = tiny_data["theta"]
        q_vec = tiny_data["cap_rates"]

        from scipy.special import expit

        p_bg = expit(X @ theta)
        I = infor_mat(X, p_bg, q_vec)

        eigenvalues = np.linalg.eigvalsh(I)
        assert np.all(eigenvalues > 0)


class TestLossGradPen:
    def test_returns_array(self, tiny_data):
        X = tiny_data["X"]
        theta = tiny_data["theta"]
        q_vec = tiny_data["cap_rates"]

        from scipy.special import expit

        p_bg = expit(X @ theta)
        wii_sqrt_X = _compute_wii_sqrt_X(X, p_bg, q_vec)
        inf_m = wii_sqrt_X.T @ wii_sqrt_X

        pen = loss_grad_pen(X, p_bg, q_vec, inf_m, wii_sqrt_X)
        assert pen is not None
        assert pen.shape == (X.shape[1],)


class TestLossFunStar:
    def test_includes_firth_penalty(self, tiny_data):
        X = tiny_data["X"]
        theta = tiny_data["theta"]
        q_vec = tiny_data["cap_rates"]
        y_vec = tiny_data["pic_matrix"][0]

        from scipy.special import expit

        p_bg = expit(X @ theta)

        base = loss_fun(p_bg, q_vec, y_vec)
        star = loss_fun_star(X, p_bg, q_vec, y_vec)
        # The regularized loss includes a Firth penalty term (0.5 * log|I|)
        # which can be positive or negative depending on the information matrix
        assert np.isfinite(star)
        assert star != base  # Should differ from base loss


class TestIRLS:
    def test_returns_valid_result(self):
        """Test IRLS returns valid parameters and a convergence status."""
        rng = np.random.RandomState(42)
        n_cells = 30
        X = np.column_stack([np.ones(n_cells), rng.randn(n_cells)])
        q_vec = rng.uniform(0.3, 0.8, n_cells)
        y_vec = rng.binomial(1, 0.4, n_cells).astype(np.float64)
        theta_init = np.array([0.05, 0.05])

        result = irls_iter(y_vec, X, theta_init, q_vec)
        assert len(result) == X.shape[1] + 1

        conv_status = int(result[-1])
        theta_est = result[:-1]
        # Should return a valid convergence status
        assert conv_status in [s.value for s in ConvergenceStatus]
        assert np.all(np.isfinite(theta_est))

    def test_converges_with_good_data(self):
        """Test IRLS converges with well-conditioned data."""
        rng = np.random.RandomState(100)
        n_cells = 100
        X = np.column_stack([np.ones(n_cells), rng.randn(n_cells) * 0.5])
        q_vec = np.full(n_cells, 0.5)
        # Generate data from the model
        true_theta = np.array([-0.5, 0.3])
        p_true = 1.0 / (1.0 + np.exp(-(X @ true_theta)))
        y_vec = rng.binomial(1, p_true * q_vec).astype(np.float64)
        theta_init = np.array([0.05, 0.05])

        result = irls_iter(y_vec, X, theta_init, q_vec)
        conv_status = int(result[-1])
        assert conv_status == ConvergenceStatus.CONVERGED

    def test_null_holds_zero(self, tiny_data):
        X = tiny_data["X"]
        q_vec = tiny_data["cap_rates"]
        y_vec = tiny_data["pic_matrix"][0]
        theta_init = np.array([0.05, 0.0])
        hold_zero = np.array([1])

        result = irls_iter_null(y_vec, X, theta_init, hold_zero, q_vec)
        theta_est = result[:-1]
        # Parameter at index 1 should remain zero
        assert theta_est[1] == 0.0

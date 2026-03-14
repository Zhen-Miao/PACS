"""Tests for the Cauchy Combination Test module."""

import numpy as np

from pacs.cct import cauchy_combination_test


class TestCauchyCombinationTest:
    def test_basic(self):
        p_values = np.array([[0.05, 0.1, 0.5], [0.9, 0.8, 0.7]])
        result = cauchy_combination_test(p_values)

        assert result.shape == (2,)
        assert np.all((result >= 0) & (result <= 1))

    def test_small_pvalues(self):
        """Very small p-values should produce very small combined p-values."""
        p_values = np.array([[1e-20, 1e-18, 1e-15]])
        result = cauchy_combination_test(p_values)

        assert result.shape == (1,)
        assert result[0] < 1e-10

    def test_large_pvalues(self):
        """All large p-values should produce a large combined p-value."""
        p_values = np.array([[0.9, 0.8, 0.95]])
        result = cauchy_combination_test(p_values)

        assert result[0] > 0.5

    def test_single_test(self):
        """With a single test, combined p-value should equal the input."""
        p_values = np.array([[0.05]])
        result = cauchy_combination_test(p_values)
        np.testing.assert_allclose(result[0], 0.05, atol=1e-10)

    def test_mixed_pvalues(self):
        """One very small p-value should dominate."""
        p_values = np.array([[1e-20, 0.5, 0.9]])
        result = cauchy_combination_test(p_values)
        assert result[0] < 0.01

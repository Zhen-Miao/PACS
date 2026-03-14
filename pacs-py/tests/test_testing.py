"""Tests for the testing module (main PACS test functions)."""

import numpy as np
import pytest

from pacs.testing import pacs_test, pacs_test_logit


class TestPacsTestLogit:
    def test_returns_result(self, small_data):
        result = pacs_test_logit(
            metadata=small_data["metadata"],
            formula_full="~ cell_type",
            formula_null="~ 1",
            pic_matrix=small_data["pic_matrix"],
            cap_rates=small_data["cap_rates"],
        )

        assert result.p_values.shape == (small_data["n_features"],)
        assert np.all(
            (result.p_values >= 0) & (result.p_values <= 1)
            | np.isnan(result.p_values)
        )
        assert result.convergence.shape == (small_data["n_features"], 2)

    def test_dars_have_smaller_pvalues(self, small_data):
        result = pacs_test_logit(
            metadata=small_data["metadata"],
            formula_full="~ cell_type",
            formula_null="~ 1",
            pic_matrix=small_data["pic_matrix"],
            cap_rates=small_data["cap_rates"],
        )

        n_dars = small_data["n_dars"]
        dar_pvals = result.p_values.values[:n_dars]
        nondar_pvals = result.p_values.values[n_dars:]

        # On average, DAR p-values should be smaller
        valid_dar = dar_pvals[np.isfinite(dar_pvals)]
        valid_nondar = nondar_pvals[np.isfinite(nondar_pvals)]

        if len(valid_dar) > 0 and len(valid_nondar) > 0:
            assert np.median(valid_dar) < np.median(valid_nondar)


class TestPacsTest:
    def test_auto_selection(self, small_data):
        result = pacs_test(
            pic_matrix=small_data["pic_matrix"],
            metadata=small_data["metadata"],
            formula_full="~ cell_type",
            formula_null="~ 1",
            cap_rates=small_data["cap_rates"],
            verbose=False,
        )

        assert result.p_values.shape == (small_data["n_features"],)
        assert result.convergence.shape == (small_data["n_features"], 2)

    def test_input_validation(self, small_data):
        # Mismatched dimensions
        with pytest.raises(ValueError, match="cap_rates"):
            pacs_test(
                pic_matrix=small_data["pic_matrix"],
                metadata=small_data["metadata"],
                formula_full="~ cell_type",
                formula_null="~ 1",
                cap_rates=small_data["cap_rates"][:5],
                verbose=False,
            )

    def test_sparse_input(self, small_data):
        import scipy.sparse as sp

        pic_sparse = sp.csr_matrix(small_data["pic_matrix"])
        result = pacs_test(
            pic_matrix=pic_sparse,
            metadata=small_data["metadata"],
            formula_full="~ cell_type",
            formula_null="~ 1",
            cap_rates=small_data["cap_rates"],
            verbose=False,
        )

        assert result.p_values.shape == (small_data["n_features"],)

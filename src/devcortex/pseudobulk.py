"""Pseudobulk aggregation of a (cells x features) sparse count matrix by cell state."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def state_labels(obs: pd.DataFrame, cell_type_col: str, stage_col: str) -> pd.Series:
    """'cell_type|stage' label for each nucleus."""
    return obs[cell_type_col].astype(str) + "|" + obs[stage_col].astype(str)


def aggregate_pseudobulk(x, labels: pd.Series, keep_cells: pd.Index | None = None) -> pd.DataFrame:
    """Sum counts per label. ``x`` is (n_cells x n_features), sparse or dense.

    Returns a features x states DataFrame of raw counts (float64). Rows follow the
    column order of ``x``; the caller attaches feature ids.
    """
    labels = labels.reset_index(drop=True) if keep_cells is None else labels.loc[keep_cells]
    if keep_cells is not None:
        x = x[np.asarray(keep_cells)]
    cats = pd.Categorical(labels)
    n = len(cats)
    ind = sp.csr_matrix((np.ones(n), (np.arange(n), cats.codes)), shape=(n, len(cats.categories)))
    xs = x.tocsr() if sp.issparse(x) else sp.csr_matrix(x)
    agg = (ind.T @ xs).T          # features x states
    return pd.DataFrame(agg.toarray(), columns=list(cats.categories))

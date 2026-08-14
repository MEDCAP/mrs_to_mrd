"""
Truncated-SVD (low-rank) denoising of complex MRS signal matrices.

Following Francischello et al., NMR Biomed. 2021;34:e4285.

The columns of the matrix are the complex time-domain FIDs of a time series of
spectra. In the absence of noise a series of metabolite spectra with fixed line
shapes has rank ~= number of metabolites; additive Gaussian noise makes it full
rank. Keeping only the first `rank` singular values recovers the low-rank signal
subspace and discards the noise subspace.

Denoising operates on the raw complex signal, so no phase correction is required
beforehand. Pure numpy - no mrd dependency, so this is usable and testable on its
own.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DenoiseResult:
    matrix: np.ndarray           # (m, n) complex, denoised
    singular_values: np.ndarray  # (min(m, n),) float64, of the *original* matrix
    rank: int                    # number of singular values retained


def gavish_donoho_rank(singular_values: np.ndarray, shape: tuple) -> int:
    """
    Gavish-Donoho optimal hard threshold for a matrix with unknown noise level.

    thresh = omega(beta) * median(S), with beta = min(m, n) / max(m, n).
    Args:
        - singular_values: the singular values of the matrix, descending
        - shape: (m, n) of the matrix they came from
    Returns:
        - the number of singular values above the threshold, at least 1
    """
    m, n = shape
    beta = min(m, n) / max(m, n)
    omega = 0.56 * beta ** 3 - 0.95 * beta ** 2 + 1.82 * beta + 1.43
    thresh = omega * np.median(singular_values)
    return max(1, int(np.count_nonzero(singular_values > thresh)))


def denoise_svd(matrix: np.ndarray, rank: int = None) -> DenoiseResult:
    """
    Low-rank approximation M_hat = U @ diag(S_r) @ Vh of a complex signal matrix.

    Args:
        - matrix: complex ndarray of shape (nsamples, nspectra); each column is a FID
        - rank: number of singular values to retain. If None, estimated with
          `gavish_donoho_rank`.
    Returns:
        - DenoiseResult(matrix=denoised, singular_values=S, rank=rank used)
    """
    U, S, Vh = np.linalg.svd(matrix, full_matrices=False)
    if rank is None:
        rank = gavish_donoho_rank(S, matrix.shape)
    rank = min(rank, len(S))
    S_trunc = S.copy()
    S_trunc[rank:] = 0
    return DenoiseResult(matrix=(U * S_trunc) @ Vh,
                         singular_values=S.astype(np.float64),
                         rank=int(rank))

import numpy as np
from scipy.linalg import eigh

def laplacian_eigendecomposition(A, k=None):
    """
    Perform Laplacian eigen-decomposition on adjacency matrix A.
    
    Parameters:
    - A (ndarray): Adjacency matrix of shape (N, N)
    - k (int or None): Number of smallest eigenvectors to return (optional)

    Returns:
    - eigenvalues (ndarray): Array of eigenvalues (sorted ascending)
    - eigenvectors (ndarray): Matrix of eigenvectors (each column is an eigenvector)
    """
    # Degree matrix
    D = np.diag(np.sum(A, axis=1))

    # Unnormalized Laplacian
    L = D - A

    # Compute eigenvalues and eigenvectors
    # eigh is used for symmetric (Hermitian) matrices like L
    eigenvalues, eigenvectors = eigh(L)

    if k is not None:
        return eigenvalues[:k], eigenvectors[:, :k]
    return eigenvalues, eigenvectors

# Example usage
A = np.array([
    [0, 1, 1, 0],
    [1, 0, 1, 0],
    [1, 1, 0, 0],
    [0, 0, 0, 0]
])

eigvals, eigvecs = laplacian_eigendecomposition(A)
print("Eigenvalues:", eigvals)
print("Eigenvectors:\n", eigvecs)

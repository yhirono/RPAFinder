import numpy as np
import sympy as sp
from scipy.linalg import null_space


def proj_mat(xsub, xall):
    mask = np.isin(xall, xsub)
    return np.diag(mask.astype(int))


def proj_mat_comp(xsub, xall):
    I = np.identity(len(xall))
    return I - proj_mat(xsub, xall)


def dimker(mat):
    nullspace = null_space(np.array(mat))
    return nullspace.shape[1]


def dimcoker(mat):
    nullspace = null_space(np.array(mat).T)
    return nullspace.shape[1]


def get_intersection(l1, l2):
    """
    Calculates the intersection of two lists of vectors.

    Args:
        l1 (list): The first list of vectors.
        l2 (list): The second list of vectors.

    Returns:
        numpy.ndarray: The intersection of the two lists of vectors.
    """
    n_l1, n_l2 = len(l1), len(l2)
    if n_l1 == 0 or n_l2 == 0:
        return np.array([])
    A = np.array(l1).transpose()
    B = np.array(l2).transpose()
    C = np.hstack((A,-B))
    U = null_space(C)
    return (A @ U[:n_l1]).transpose()


def get_intersection_symb(l1, l2):
    """
    Calculates the intersection of two symbolic matrices.

    Parameters:
    l1 (numpy.ndarray): The first symbolic matrix.
    l2 (numpy.ndarray): The second symbolic matrix.

    Returns:
    numpy.ndarray: The intersection of the two symbolic matrices.
    """
    n_l1, n_l2 = l1.shape[0], l2.shape[0]

    if n_l1 == 0 and n_l2 == 0:
        return sp.zeros(0,0)
    elif n_l1 == 0 or n_l1==0:
        return sp.zeros(0, l1.shape[1])

    A = sp.Matrix(l1).T
    B = sp.Matrix(l2).T
    C = A.row_join(-B)
    U = C.nullspace()
    if not U:
        return sp.zeros(0, l1.shape[1])

    U_matrix = sp.Matrix.hstack(*U)
    intersection = A * U_matrix[:n_l1, :]
    return intersection.T


def get_orthogonal_complement(l1):
    """
    Calculates the orthogonal complement of a given vectors.

    Parameters:
    l1 (list or numpy.ndarray): The input vectors

    Returns:
    numpy.ndarray: The orthogonal complement of the input vectors.
    """
    A = np.array(l1)
    if A.shape[0] == 0:
        return np.identity(A.shape[1])
    U = null_space(A)
    return U.transpose()


def get_orthogonal_complement_symb(l1):
    """
    Compute the orthogonal complement of a given list of vectors.

    Parameters:
    l1 (list or matrix): The input vectors.

    Returns:
    matrix: The orthogonal complement of the input vectors.
    """
    A = sp.Matrix(l1)
    if A.rows == 0:
        return sp.eye(A.cols)
    U = A.nullspace()
    if not U:
        return sp.zeros(0, A.cols)
    U_matrix = sp.Matrix.hstack(*U)
    return U_matrix.T

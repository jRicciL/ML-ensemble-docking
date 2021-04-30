import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpltcols
import matplotlib.patches as mpatches
from typing import Tuple

def cMDS(D: np.ndarray, 
         is_similarity: bool = False
         ) -> Tuple:

    '''
    Computes Classical Multidimensional Scaling from a given distance, or 
    similarity, matrix D.

    Parameters
    ----------
    D : np.ndarray
        A distance, or similarity, matrix (squared matrix)
    is_similarity: bool
        Determines if D is a similarity matrix (True) or a
        distance matrix (False)

    Returns
    -------
    Tuple
        F.T: np.ndarray
            transposed configuration matrix F
        D_sq: np.ndarray
            squared distance matrix D
        B: np.ndarray
            double centering matrix = -0.5*J*D^2*J  
        e_vals: np.array
            eigenvalues of B

    Modified from: http://www.nervouscomputer.com/hfs/cmdscale-in-python/
    '''

    assert D.shape[0] == D.shape[1]

    if is_similarity:
        D = 1 - D # If D is a similarity matrix, convert it to a distance matrix
    # Number of samples                                                             
    n = len(D)
    # Create the Squared proximity matrix
    D_sq = D**2
    # Generate the Centering matrix: J
    # Defined as the identity matrix I_n minums an nxn all-ones matrix
    e1 = np.ones((n,1))
    m = ( np.ones((n, 1)) / n)
    I = np.eye(n) # Identity matrix I_n
    J = I - e1.dot(np.transpose(m)) # Centering matrix J
    # Generate B = -0.5*J*D^2*J                                                                                 
    B = -(1/2) * J.dot(D_sq).dot(J)

    # Compute the eigendecomposition of B                                                                         
    e_vals, e_vecs = np.linalg.eigh(B)

    # Sort eigenvalues (descending)
    idx   = np.argsort(e_vals)[::-1] # Get the sorting indices
    e_vals = e_vals[idx] # Sort eigenvectors
    e_vecs = e_vecs[:,idx]
    w, = np.where(e_vals > 0) # To keep only positive eigenvalues
    e_vals = e_vals[w]
    # Configuration matrix F = M^(-1/2)*E_m*L_m^(1/2)  
    M_neg_sq = np.diag(m[:,0]**(-0.5)) # M_neg_sq is the  M^(-1/2) matrix
    L_m_sq  = np.diag(np.sqrt(e_vals[w])) #  Lambda_sq (L_m^(1/2)) ia a diagonal matrix with the m positive eigenvalues
    E_m = e_vecs[:,w] # E_m  = sorted eigenvectors (positives only)
    # F  = M_neg_sq.dot(E_m.dot((L_m_sq)))  # F = M^(-1/2)*E_m*L_m^(1/2)
    F = E_m.dot((L_m_sq))
    # Return the computed matrices
    return F.T, D_sq, B, e_vals

def cMDS_proj(cMDS_obj: Tuple, 
              sup_point: np.ndarray) -> np.ndarray:
    '''
    Computes Classical Multidimensional Scaling projection given
    a cMDS_obj (reference space) and a matrix of distances between
    supplementary and reference points.

    Parameters
    ----------
    cMDS_ob : Tuple
        A tuple with `F.T, D_sq, B, e_vals` values, returned from the
        cMDS function
    sup_point: np.ndarray
        A distance matrix of size nxm, where n is the number of 
        reference samples and m is the number of supplementary samples 

    Returns
    -------
    F_proj: np.ndarray
        A matrix with supplementary samples projected over the reference space

    Modified from: http://www.nervouscomputer.com/hfs/cmdscale-in-python/
    '''
    F = cMDS_obj[0].T # Matriz F de scores obtenida con cMDS a partir de los valores "Activos"
    D_sq = cMDS_obj[1] # Matrix Delta^2
    e_vals = cMDS_obj[3] # Eigenvalores positivos obtenidos del B de los valores "Activos"
    n = len(D_sq)
    e1 = np.ones((n, 1))
    m = len(sup_point)
    e1_sup = np.ones((m, 1))
    m = ( np.ones((n, 1)) / n)
    I = np.eye(n) # Matriz de identidad I_n
    J = I - e1.dot(np.transpose(m)) # Mariz de centrado
    # MDS out-of-smaple: Adición de un punto suplementario
    a_out = sup_point**2
    # D_sup_2 es un vector con los valores de RMSD del nuevo punto frente a los "Activos"
    # Se calcula B_sup como -0.5*J*(a_out - D^2*(m1^t))
    B_sup = -0.5 * J.dot(a_out.T - D_sq.dot(m.dot(np.transpose(e1_sup))))
    # Finalmente se calcula F_sup = (B_sup^t)*F*L_inv
    # Donde L_inv es la matriz diagonal con la inversa de los eingenvalores
    L_m_inv = np.diag( 1 / e_vals )
    F_proj = np.transpose(B_sup).dot( F ).dot( L_m_inv)
    return F_proj



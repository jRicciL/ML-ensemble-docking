from statsmodels.stats.libqsturng import psturng, qsturng
from scipy import stats
import itertools as it
import pandas as pd
import numpy as np

def _col_sig_p_values(val):
    color = 'red' if val < 0.05 else 'none'
    return 'color: %s' % color

def _col_H0_rejected(val):
    color = 'red' if val else 'green'
    return 'color: %s' % color

def get_rankings(df):
    return df.rank(axis=1, ascending=False)

def get_R(df):
    rankings = get_rankings(df)
    R_vec = rankings.mean(axis=0)
    return R_vec

def get_R2_sum(df):
    return get_R(df).apply(lambda x: x**2).sum()

def get_chi_friedman(df):
    N, k = df.shape
    first = (12 * N) / (k*(k + 1))
    second = get_R2_sum(df) - ( k*((k+1)**2) / 4 )
    return first * second

def friedmanTest(df):
    chi = get_chi_friedman(df)
    deg_f = df.shape[1] - 1
    p_value = stats.chi2.sf(chi, df=deg_f)
    
    res = pd.DataFrame({'chi^2': chi, 'dof': deg_f, 'p_value': p_value}, index=['Friedman Test'])
    return res


def friedman_imanDavenportTest(df):
    N, k = df.shape
    deg_f = df.shape[1] - 1
    chi = get_chi_friedman(df)
    F = ((N -1) * chi) / ( (N* (k - 1)) - chi )
    deg_f_1 = k -1
    deg_f_2 = (k-1)*(N-1)
    p_value= stats.f.sf(F, deg_f_1, deg_f_2)
    
    res = pd.DataFrame({'F': F, 'dof1': deg_f_1, 'dof2': deg_f_2, 'p_value': p_value}, index=['Friedman Test'])
    return res


# Pos hoc analysis
def get_critical_values_nemenyi(alpha, k, N=np.inf):
    q = qsturng(1 - alpha, k, N) / np.sqrt(2)
    return q

def get_cd_nemeyi(alpha, k, N = np.inf):
    q = get_critical_values_nemenyi(alpha, k, N)
    cd = q * np.sqrt(k*(k+1) / (6*N))
    return cd

def get_rank_diff(R_i, R_j):
    diff = abs(R_i - R_j)
    return diff

def get_q(R_i, R_j, k, N):
    diff = get_rank_diff(R_i, R_j)
    q = diff /  np.sqrt( k*(k+1) / (6*N) )
    return q

def get_p_value(q, k, N=np.inf):
    z = q * np.sqrt(2)
    p = psturng(z, k, N)
    return p

def get_data_frame(narray, names, fill=1):
    narray = narray + narray.T
    np.fill_diagonal(narray, fill)
    df = pd.DataFrame(narray, index=names, columns=names)\
                        .reindex(sorted(names, reverse=False), axis=0)\
                        .reindex(sorted(names, reverse=False), axis=1)
    return df
    

def pairwise_nemenyi(df, alpha=0.05):
    N, k = df.shape
    p_vals = np.zeros((k, k))
    diff_vals = np.zeros((k, k))
    combinations = it.combinations(range(k), 2)
    
    tri_lower = np.tril_indices(k, -1)
    tri_upper = np.triu_indices(k, 1)
    
    # Ranks
    R = get_R(df)
    # Critical difference CD
    CD = get_cd_nemeyi(alpha, k, N)
    
    for i, j in combinations:
        diff_vals[i, j] = get_rank_diff(R[i], R[j])
        p_vals[i, j] = get_p_value( get_q(R[i], R[j], k, N), k )
    
    names = df.columns
    p_vals = get_data_frame(p_vals, names)
    signif = diff_vals > CD
    signif = get_data_frame(signif, names, fill=False)
    
    return p_vals, signif
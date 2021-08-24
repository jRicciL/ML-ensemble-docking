import pandas as pd
import numpy as np
from typing import Iterable

def get_mean_score(df: pd.DataFrame, 
                   lower_is_best: bool = True):
    '''Get a dataframe of m*n values and 
       returns a numpy array of means per row'''
    df_means = df.mean(axis = 1).to_numpy()
    if lower_is_best:
        df_means *= -1
        # times -1 because roc_auc_score assumes higher is better
    return df_means


def get_min_score(df: pd.DataFrame, 
                  lower_is_best: bool = True):
    '''Get a dataframe of m*n values and 
       returns a numpy array of min scores per row'''
    df_best = df.min(axis = 1).to_numpy()
    if lower_is_best:
        df_best *= -1
        # times -1 because roc_auc_score assumes higher is better
    return df_best


def geo_mean_overflow(iterable: Iterable):
    a = np.log(iterable)
    return np.exp(a.sum()/len(a))


def get_geom_mean_score(df: pd.DataFrame, 
                        lower_is_best: bool = True):
    '''Get a dataframe of m*n values and 
       returns a numpy array of euclidean norm scores per row'''
    if lower_is_best:
        df = df.copy() * -1
    df_geom_mean = df.apply(geo_mean_overflow, axis=1)
    return df_geom_mean.to_numpy()
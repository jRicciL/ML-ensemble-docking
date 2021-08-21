# Nested vs non-nested cross-validation
# Based on: https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html

# Filename with the results
file_name = '../4_Ensemble_docking_results/df_DkSc_results_COCRYS_DEKOIS_DUD.pkl'

import pickle
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import \
        GridSearchCV, RandomizedSearchCV, \
        cross_val_score, StratifiedKFold
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings('ignore')

# Random trials
PROT_NAME  = 'cdk2'
MODEL      = 'XGB'
N_CPUS     = 16
NUM_TRIALS = 30
k_CV_INNER = 4
k_CV_OUTER = 4
ESTIMATOR  = XGBClassifier 
OPTIMIZER  = RandomizedSearchCV 

# Import the data
df_dk_res = pd.read_pickle(file_name)

# Extract the features columns: Docking scores
X = df_dk_res.drop('activity', axis = 1).values
y = df_dk_res['activity'].values

# Previously optimized params 
# CDK2 - XGBoost 
# Estimator with pre-optimized hyprms
pre_optimized_hyparams = {
    'subsample'       : 0.5, 
    'n_estimators'    : 200, 
    'max_depth'       : 10, 
    'learning_rate'   : 0.1,
    'alpha'           : 0.5,
    'gamma'           : 1, 
    'colsample_bytree': 1,
}

# GRID HYPERPARAMETERS
grid_hyprms = {
   'n_estimators'    : [200, 300, 500],
   'max_depth'       : [3, 5, 10, 20],
   'gamma'           : [0.01, 0.1, 0.5, 1],
   'learning_rate'   : [0.05, 0.1],
   'subsample'       : [0.3, 0.5, 0.6],
   'alpha'           : [0.01, 0.1, 0.5, 1],
   'colsample_bytree': [0.3, 0.5, 1]
}

# Arrays to store the scores per repetition (means)
default_hypms_scores = np.zeros(NUM_TRIALS)
preOpt_hypms_scores  = np.zeros(NUM_TRIALS)
nested_scores        = np.zeros(NUM_TRIALS)
non_nested_scores    = np.zeros(NUM_TRIALS)


# Arrays to stroe the scores per cv validation test inside repetition
ind_default_hypms_scores = np.zeros(NUM_TRIALS*k_CV_OUTER)
ind_preOpt_hypms_scores  = np.zeros(NUM_TRIALS*k_CV_OUTER)
ind_nested_scores        = np.zeros(NUM_TRIALS*k_CV_OUTER)


# ********* Loop for the trial *********
for i in range(NUM_TRIALS):
    # Here different cross-validation techniques could be applied
   
    # ********** Perform the splittings **********
    # Outer splitting
    outer_cv = StratifiedKFold(n_splits = k_CV_OUTER, 
                               shuffle  = True,
                               random_state = i)
    # Inner splitting
    inner_cv = StratifiedKFold(n_splits = k_CV_INNER,
                               shuffle  = True,
                               random_state = i)

    
    # ********** Estimator with default hyprms **********
    estimator_DEFAULT   = ESTIMATOR(
            eval_metric = 'logloss',
            use_label_encoder =  False,
            n_jobs  = N_CPUS
            )
    default_hypms_score = cross_val_score(
        estimator = estimator_DEFAULT,
        X  = X, 
        y  = y,
        cv = outer_cv,
        scoring = 'roc_auc',
    )
    default_hypms_scores[i] = default_hypms_score.mean()
    ind_default_hypms_scores[i*k_CV_OUTER:(i+1)*k_CV_OUTER] = default_hypms_score
    
    
    # ********** Estimator with pre optimized hyprms **********
    estimator_PREOPTIMIZED = ESTIMATOR(
            eval_metric = 'logloss',
            use_label_encoder =  False,
            n_jobs  = N_CPUS,
            **pre_optimized_hyparams)
    preOpt_hypms_score = cross_val_score(
        estimator = estimator_PREOPTIMIZED,
        X  = X, 
        y  = y,
        cv = outer_cv,
        scoring = 'roc_auc',
    )
    preOpt_hypms_scores[i] = preOpt_hypms_score.mean()
    ind_preOpt_hypms_scores[i*k_CV_OUTER:(i+1)*k_CV_OUTER] = preOpt_hypms_score
    
    
    # ********** Estimator with optimized hyprms inside outer loop **********
    estimator = ESTIMATOR(
            eval_metric = 'logloss',
            use_label_encoder =  False,
            n_jobs = N_CPUS,
            )
    clf_gs = RandomizedSearchCV(
                       estimator = estimator,
                       param_distributions = grid_hyprms,
                       cv = inner_cv,
                       scoring = 'roc_auc',
                       n_iter = 50
                )
    clf_gs.fit(X, y)
    non_nested_scores[i] = clf_gs.best_score_


    # ********** Nested CV with parameter optimization **********
    # Inside each fold of the cross_val_score perform a GridSearch
    nested_score = cross_val_score(
        estimator = clf_gs,
        X  = X, 
        y  = y,
        cv = outer_cv,
        scoring = 'roc_auc',
    )
    nested_scores[i] = nested_score.mean()
    ind_nested_scores[i*k_CV_OUTER:(i+1)*k_CV_OUTER] = nested_score

    # Save the results per repetition

    # repets*k_outers dataframe
    df_1 = pd.DataFrame({
        'nested': ind_nested_scores,
        'non_nested': np.repeat(non_nested_scores, k_CV_INNER),
        'preOptHpm': ind_preOpt_hypms_scores,
        'defHpm': ind_default_hypms_scores
    })
    df_1.to_csv(f'./DF1_{PROT_NAME}_{MODEL}_{NUM_TRIALS}reps_{k_CV_OUTER}x{k_CV_INNER}nCV_results.csv')

    # repeats
    df_2 = pd.DataFrame({
        'nested': nested_scores,
        'non_nested': non_nested_scores,
        'pre_optimized_hyprms': preOpt_hypms_scores,
        'default_hyprms': default_hypms_scores
    })
    df_2.to_csv(f'./DF2_MEANS_{PROT_NAME}_{MODEL}_{NUM_TRIALS}reps_{k_CV_OUTER}x{k_CV_INNER}nCV_results.csv')

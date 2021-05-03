# first line: 1
@memory.cache
def nk_rep_cross_validation(
          estimators, X, y,
          metrics,
          n_splits=2, 
          n_repeats=5,
          random_state=None, 
          shuffle=True,  y_preds_return=False):
    # Compute the Stratified K folds
    cv = RepeatedStratifiedKFold(
                         n_splits=n_splits,
                         n_repeats=n_repeats,
                         random_state=random_state)
    splits = [*cv.split(X, y)]
    
    results, y_preds_dict = _do_replicates(splits, 
                                estimators, X, y, 
                                metrics)
    
    df_res = _format_results_to_df(metrics, results, 
                                   n=n_splits*n_repeats)
    
    if y_preds_return:
        return df_res, y_preds_dict, splits  # Test idx
    else:
        return df_res 

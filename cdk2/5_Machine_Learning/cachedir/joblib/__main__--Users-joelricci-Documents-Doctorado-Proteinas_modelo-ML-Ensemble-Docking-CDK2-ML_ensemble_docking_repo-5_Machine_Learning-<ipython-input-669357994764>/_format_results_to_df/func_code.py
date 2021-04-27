# first line: 1
@memory.cache
def _format_results_to_df(metrics, results, n):
        # Format into a dataframe
    # Create the metric names and repeat them 
    n_metrics = len(metrics)
    index_names = [*metrics.keys()]*n
    
    # convert to a dataframe
    df_res = pd.DataFrame(
        results, 
        index= pd.MultiIndex.from_tuples(
            zip(index_names,
                np.repeat(range(n), n_metrics))
        ))
    df_res = df_res.sort_index()
    
    return df_res

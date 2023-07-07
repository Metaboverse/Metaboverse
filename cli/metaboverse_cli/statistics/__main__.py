import pandas as pd 
from scipy.stats import ttest_ind 
import sys
import json
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def __main__(args_dict):

    # Read data from stdin
    input_data = json.loads(sys.stdin.read())

    # Convert lists to numpy arrays
    arr1 = np.array(input_data['array1'])
    arr2 = np.array(input_data['array2'])

    # Perform t-test
    t_stats, p_values = stats.ttest_ind(arr1, arr2, axis=1)

    # Perform Benjamini/Hochberg adjustment
    _, p_values_adjusted, _, _ = multipletests(p_values, method='fdr_bh')

    # Calculate the sample sizes, means and standard deviations
    n1 = arr1.shape[1]
    n2 = arr2.shape[1]
    mean1 = np.mean(arr1, axis=1)
    mean2 = np.mean(arr2, axis=1)
    std1 = np.std(arr1, axis=1, ddof=1)
    std2 = np.std(arr2, axis=1, ddof=1)

    # Calculate the standard error
    se1 = std1 / np.sqrt(n1)
    se2 = std2 / np.sqrt(n2)

    # Degrees of freedom
    df = n1 + n2 - 2

    # Calculate the 90%, 95% and 99% confidence intervals for both arrays
    conf_intervals_90 = stats.t.interval(0.9, df, mean1 - mean2, np.sqrt(se1**2 + se2**2))
    conf_intervals_95 = stats.t.interval(0.95, df, mean1 - mean2, np.sqrt(se1**2 + se2**2))
    conf_intervals_99 = stats.t.interval(0.99, df, mean1 - mean2, np.sqrt(se1**2 + se2**2))

    # Output adjusted p-values
    # Output p-values and confidence intervals
    print(json.dumps({
        "unadjusted": p_values.tolist(), 
        "adjusted": p_values_adjusted.tolist(),
        "conf_intervals_90": np.array(conf_intervals_90).tolist(),
        "conf_intervals_95": np.array(conf_intervals_95).tolist(),
        "conf_intervals_99": np.array(conf_intervals_99).tolist()
    }))






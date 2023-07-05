import pandas as pd 
import ttest_ind from scipy.stats
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
    _, p_values = stats.ttest_ind(arr1, arr2, axis=1)

    if args_dict['type'] == 'ttest':
        # Output p-values
        print(json.dumps(p_values.tolist()))
    else:
        # Perform Benjamini/Hochberg adjustment
        _, p_values_adjusted, _, _ = multipletests(p_values, method='fdr_bh')

        # Output adjusted p-values
        print(json.dumps(p_values_adjusted.tolist()))






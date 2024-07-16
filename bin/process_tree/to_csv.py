"""
Prep inputs for tree processing.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


##


# Args
path_input = sys.argv[1]

# Main
def main():

    # AD and DP
    AD = load_npz(os.path.join(path_input, 'AD.npz')).astype(np.int16)
    DP = load_npz(os.path.join(path_input, 'DP.npz')).astype(np.int16)
    meta = pd.read_csv(os.path.join(path_input, 'meta.csv'), index_col=0)
    var_meta = pd.read_csv(os.path.join(path_input, 'var_meta.csv'), index_col=0)
    pd.DataFrame(AD.A, index=meta.index, columns=var_meta.index).to_csv('AD.csv')
    pd.DataFrame(DP.A, index=meta.index, columns=var_meta.index).to_csv('DP.csv')

# Run 
if __name__ == '__main__':
    main()


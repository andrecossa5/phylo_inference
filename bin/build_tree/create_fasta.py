#!/usr/bin/python

# create_fasta script

import os
import sys
from mito_utils.utils import *
from mito_utils.preprocessing import *


# Main
def main():

    # os.chdir('/Users/IEO5505/Desktop/MI_TO/phylo_inference/scratch')
    # path_afm = 'afm_new.h5ad'

    path_afm = sys.argv[1]
    afm = sc.read(path_afm)
    d = AFM_to_seqs(afm)
    with open('genotypes.fa', 'w') as f:
        for cell in d:
            f.write(f'>{cell}\n{d[cell]}\n')


if __name__ == '__main__':
    main()

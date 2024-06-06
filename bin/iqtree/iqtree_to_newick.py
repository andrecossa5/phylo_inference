#!/usr/bin/python

"""
Convert a .iqtree object to a newick string.
"""

import os


##


def load_newick_from_iqtree():
    newick_format = ""
    with open('sequences.fasta.iqtree', 'r') as file:
        for line in file:
            if line.strip().startswith("("):  # Newick trees start with '('
                newick_format = line.strip()
                break
    return newick_format


##

if __name__ == '__main__':

    # Read 
    newick_format_string = load_newick_from_iqtree()
    # Write
    with open(f'tree.newick', 'w') as file:
        file.write(newick_format_string)


##
#!/usr/bin/python

"""
Convert a .iqtree object to a newick string.
"""


##


def main():

    newick_format = ""
    with open('sequences.fasta.iqtree', 'r') as file:
        for line in file:
            if line.strip().startswith("("):  # Newick trees start with '('
                newick_format = line.strip()
                break

    with open(f'tree.newick', 'w') as file:
        file.write(newick_format)


##


if __name__ == '__main__':
    main()


##
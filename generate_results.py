###########################################################################
#   Generates the figures from the expression table
###########################################################################

import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import os
from os.path import join
import sys
from optparse import OptionParser
from collections import defaultdict
import numpy as np

from onto_lib_py3 import load_ontology

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    #parser.add_option("-a", "--a_descrip", action="store_true", help="This is a flat")
    parser.add_option("-o", "--out_dir", help="Directory in which to write the ouptut")
    (options, args) = parser.parse_args()

    cell_line_f = args[0]
    expression_f = args[1]
    out_dir = options.out_dir

    og = load_ontology.load('17')[0]


    cell_line_df = pd.read_csv(cell_line_f, sep='\t', index_col=0)
    expression_df = pd.read_csv(expression_f, sep='\t', index_col=0)

    ace2_expr_df = expression_df.loc[expression_df['ENSG00000130234'] > 1.0]
    tmprss2_expr_df = expression_df.loc[expression_df['ENSG00000184012'] > 1.0]    

    cell_line_df['Count'] = np.full(len(cell_line_df), 1)
    ace2_expr_df['Count'] = np.full(len(ace2_expr_df), 1) 
    tmprss2_expr_df['Count'] = np.full(len(tmprss2_expr_df), 1)

    all_cell_line_counts_df = cell_line_df.groupby(by='Cellosaurus cell line').count()
    ace2_cell_line_counts_df = ace2_expr_df.groupby(by='Cellosaurus cell line').count()
    tmprss2_cell_line_counts_df = tmprss2_expr_df.groupby(by='Cellosaurus cell line').count()
   
 
    ace2_df = ace2_cell_line_counts_df.join(all_cell_line_counts_df, lsuffix='_ACE2', rsuffix='_all')[['Count_ACE2', 'Count_all']]
    tmprss2_df = tmprss2_cell_line_counts_df.join(all_cell_line_counts_df, lsuffix='_TMPRSS2', rsuffix='_all')[['Count_TMPRSS2', 'Count_all']]
    ace2_df['Fraction'] = ace2_df['Count_ACE2'] / ace2_df['Count_all']
    tmprss2_df['Fraction'] = tmprss2_df['Count_TMPRSS2'] / tmprss2_df['Count_all']
    ace2_df['Cell line name'] = [
        og.id_to_term[x].name
        for x in ace2_df.index
    ] 
    tmprss2_df['Cell line name'] = [
        og.id_to_term[x].name
        for x in tmprss2_df.index
    ]
    ace2_df = ace2_df.sort_values(by='Fraction', ascending=False)
    tmprss2_df = tmprss2_df.sort_values(by='Fraction', ascending=False)
    print(ace2_df)
    print(tmprss2_df)
    ace2_df.to_csv(join(out_dir, 'ACE2.tsv'), sep='\t')
    tmprss2_df.to_csv(join(out_dir, 'TMPRSS2.tsv'), sep='\t')
    return    
    

if __name__ == "__main__":
    main()

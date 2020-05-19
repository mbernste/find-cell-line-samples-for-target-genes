from optparse import OptionParser
import json
from collections import defaultdict
import sys
import pandas as pd

import kallisto_quantified_data_manager_hdf5_py3 as kqdm

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    experiment_cell_line_f = args[0]
    expressed_exps_f = args[1]
    target_genes = args[2].split(',')
    out_f = options.out_file
    with open(expressed_exps_f, 'r') as f:
        expr_exps = [l.strip() for l in f]


    expr_exps, data_matrix, gene_names = kqdm.get_gene_tpms_for_experiments(expr_exps)

    cell_line_df = pd.read_csv(experiment_cell_line_f, sep='\t', index_col=0)
    expr_df = pd.DataFrame(
        data=data_matrix,
        index=expr_exps,
        columns=gene_names
    )
    expr_df = expr_df[target_genes]

    df = expr_df.join(cell_line_df)
    
    df.to_csv(out_f, sep='\t')


if __name__ == "__main__":
    main()

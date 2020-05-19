from optparse import OptionParser
import sys
import pandas as pd

import kallisto_quantified_data_manager_hdf5_py3 as kqdm

CHUNK_SIZE = 500

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()
    
    cell_line_exps_f = args[0]
    target_genes = args[1].split(',')
    out_f = options.out_file

    exps = []
    df = pd.read_csv(cell_line_exps_f, sep='\t', index_col=0)
    exps = df.index
   
    exps_in_db = list(kqdm.filter_for_experiments_in_db(exps))
    print("Found quantified results for {}/{} experiments.".format(len(exps_in_db), len(exps))) 

    genes_in_db = set(kqdm.get_all_gene_names_in_hg38_v24_kallisto())
    print(target_genes)
    assert set(target_genes) < set(genes_in_db)
    remove_genes = genes_in_db - set(target_genes)

    has_target_gene_expressed = []
    n_processed = 0
    for exp_chunk in chunks(exps_in_db, CHUNK_SIZE):
        returned_exps, data_matrix, gene_names = kqdm.get_gene_tpms_for_experiments(
           exp_chunk,
            remove_genes = remove_genes
        )
        has_target_gene_expressed = []
        for exp, vec in zip(returned_exps, data_matrix):
            if max(vec) > 1.0:
                has_target_gene_expressed.append(exp)
        n_processed += len(exp_chunk)
        print("Processed {}/{} experiments...".format(n_processed, len(exps_in_db))) 
    with open(out_f, 'w') as f:
        f.write('\n'.join(has_target_gene_expressed))

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n+1):
        yield l[i:i + n + 1] 

if __name__ == "__main__":
    main()

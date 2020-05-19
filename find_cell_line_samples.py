from optparse import OptionParser
import os
from os.path import join
from optparse import OptionParser
import subprocess
import sqlite3
from collections import defaultdict
import pandas as pd

from onto_lib_py3 import ontology_graph

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="File in which to write output")
    (options, args) = parser.parse_args()

    metasra_f = args[0]
    sra_f = args[1]
    out_f = options.out_file

    sample_to_exps = map_sample_to_experiments(sra_f)
    sample_to_cell_line = map_samples_to_cell_line(metasra_f)

    da = []
    for sample, cell_line in sample_to_cell_line.items():
        for exp in sample_to_exps[sample]:
            da.append((sample, exp, cell_line))
    df = pd.DataFrame(
        data=da,
        columns=[
            'Sample accession', 
            'Experiment accession',     
            'Cellosaurus cell line'
        ]
    )
    df = df.set_index('Experiment accession')
    df.to_csv(out_f, sep='\t')


def map_sample_to_experiments(sra_f):
    """ 
    Map each Biosample entry to its SRA experiments.
    """
    experiment_sql = """SELECT experiment_accession,
    sample_accession, study_accession FROM experiment
    """
    print("Querying database for experiment to sample mappings...")
    sample_to_exps = defaultdict(lambda: set())
    with sqlite3.connect(sra_f) as db_conn:
        c = db_conn.cursor()
        returned = c.execute(experiment_sql)
        for r in returned:
            exp = r[0]
            sample = r[1]
            study = r[2]
            sample_to_exps[sample].add(exp)
    print("done.")
    return sample_to_exps


def map_samples_to_cell_line(metasra_f):
    sql = "SELECT sample_accession, term_id FROM mapped_ontology_terms;"
    print("Querying database for sample to terms mappings...")
    sample_to_mapped_terms = defaultdict(lambda: set())
    with sqlite3.connect(metasra_f) as metasra_conn:
        metasra_c = metasra_conn.cursor()
        results = metasra_c.execute(sql)
        for r in results:
            sample = r[0]
            term_id = r[1]
            sample_to_mapped_terms[sample].add(term_id)
    sample_to_cell_lines = defaultdict(lambda: set())
    for sample, terms in sample_to_mapped_terms.items():
        for term in terms:
            if 'CVCL' in term:
                sample_to_cell_lines[sample].add(term)
    sample_to_cell_line = {}
    for sample, cell_lines in sample_to_cell_lines.items():
        if len(cell_lines) == 1:
            sample_to_cell_line[sample] = list(cell_lines)[0]
        else:
            print('Multiple samples mapped to sample {}. Excluding from results.'.format(sample))
    print("done.")
    return sample_to_cell_line



if __name__ == "__main__":
    main()

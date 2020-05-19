#

configfile: 'config.json'

rule all:
    input: 
        '{}/ACE2.tsv'.format(config['output']),
        '{}/TMPRSS2.tsv'.format(config['output'])

rule find_cell_line_samples:
    output:
        '{}/cell_line_experiments.tsv'.format(config['output'])
    run:
        cmd = 'python find_cell_line_samples.py {metasra} {sra} -o {{output}}'.format(
            metasra=config['metasra_file'],
            sra=config['sra_file']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule filter_samples_by_target_gene_expression:
    input:
        '{}/cell_line_experiments.tsv'.format(config['output'])
    output:
        '{}/experiments_target_gene_expressed.tsv'.format(config['output'])
    run:
        cmd = 'python filter_experiments_by_gene_set_max_tpm.py {{input}} {genes} -o {{output}}'.format(
            genes=','.join(config['target_genes'])
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule build_expression_table:
    input:
        all='{}/cell_line_experiments.tsv'.format(config['output']),
        filtered='{}/experiments_target_gene_expressed.tsv'.format(config['output'])
    output:
        '{}/expression_table.tsv'.format(config['output'])
    run:
        cmd = 'python build_expression_table.py {{input.all}} {{input.filtered}} {genes} -o {{output}}'.format(
            genes=','.join(config['target_genes'])
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

# Note, this final step is specific to ACE2 and TMPRSS2. Will need to generalize
# if we consider more genes.
rule generate_final_results_table:
    input:
        all='{}/cell_line_experiments.tsv'.format(config['output']),
        expr='{}/expression_table.tsv'.format(config['output'])
    output:
        '{}/ACE2.tsv'.format(config['output']),
        '{}/TMPRSS2.tsv'.format(config['output'])
    run:
        commands = [
            'python2.7 {src_dir}/generate_results.py {{input.all}} {{input.expr}} -o {out_dir}'.format(
                src_dir=SRC_DIR,
                out_dir=OUT_DIR
            )
        ]
        for c in commands:
            shell(c)



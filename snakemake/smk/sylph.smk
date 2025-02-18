SYLPH = TOOLS_LOCATIONS['sylph']['path']
SYLPH_DB = TOOLS_LOCATIONS['sylph']['db_path']
SYLPH_TAX = TOOLS_LOCATIONS['sylph']['sylph_tax']

"""
Classifies reads/assemblies using sylph
"""
rule sylph_classification:
    """
    Sketching reads and gnerate a small index
    """

    input:
        lambda wildcards: INPUT_READS
    output:
        index = ROOT / 'sylph' / 'reads.sylsp',
        TSV = ROOT / 'sylph' / 'classification.tsv'
    params:
        db = SYLPH_DB
    threads:32
    shell:
    """
    {SYLPH} sketch -r {input} -d {output} -S reads -t {threads}
    {SYLPH} profile --estimate-unknown --read-seq-id 95 --min-number-kmers 40 {params.db} {output.index} -o {output.TSV}
    """
rule sylph_clean:

    """
    clean output of sylph
    
    """
    input:
        rules.sylph_classification.TSV
    output:
        TSV = ROOT / 'sylph' / 'classification_clean.tsv'
    threads: 8
    shell:
    """
    python {snakemake.input} {snakemake.output}
    """

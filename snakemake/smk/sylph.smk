SYLPH = TOOLS_LOCATIONS['sylph']['path']
SYLPH_DB1 = TOOLS_LOCATIONS['sylph']['db_path1']
SYLPH_DB2 = TOOLS_LOCATIONS['sylph']['db_path2']
SYLPH_TAX = TOOLS_LOCATIONS['sylph']['sylph_tax']
METADATA = TOOLS_LOCATIONS['sylph']['sylph_tax_metadata']
OUTPUT_DIR= ROOT / 'sylph'

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
        TSV = ROOT / 'sylph' / 'classification.tsv'
    params:
        db1 = SYLPH_DB1,
        db2 = SYLPH_DB2,
    threads:30
    shell:
        """
        {SYLPH} profile {params.db1} {params.db2} {input}  --estimate-unknown --read-seq-id 95 --min-number-kmers 40 -t {threads} -o {output.TSV}
        """
rule sylph_taxprof:
    input:
        rules.sylph_classification.output.TSV
    output:
        TSV = ROOT / 'sylph' / 'classification_tax.tsv'
    params:
        db = METADATA,
        prefix = ROOT / 'sylph' / 'classification_tmp',
        output_dir = OUTPUT_DIR,
    threads: 8
    shell:
        """
        cd {params.output_dir}
        sylph-tax taxprof {input} -o {params.prefix} -t {params.db}
        cat classification_tmp*.sylphmpa > {output.TSV}
        """
rule sylph_clean:
    """
    generate metadata
    """
    input:
        classification_tax = rules.sylph_taxprof.output.TSV
    output:
        TSV = ROOT / 'sylph' / 'classification_clean.tsv',
        TAX_TSV = ROOT / 'sylph' / 'classification_cleaned_tax.tsv'
    threads: 2
    run:
        import pandas as pd
        import re
        species_names = []
        genus_names = []
        # abundances = []
        relative_abundances = []
        sequence_abundances = []
        with open(input.classification_tax, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if re.search(r'\|s__',line):
                    clade_name, relative_abundance, sequence_abundance = line.strip().split('\t')[:3]
                    if re.search(r'\|t__',clade_name):
                        continue
                    else:
                        species_name = clade_name.split('|')[-1].split('_')[2]
                        genus_name = clade_name.split('|')[-1].split('_')[2].split(' ')[0].split('_')[0]
                    species_names.append(species_name)
                    genus_names.append(genus_name)
                    relative_abundances.append(f'{float(relative_abundance):.5f}')
                    sequence_abundances.append(f'{float(sequence_abundance):.5f}')

        pd.Series(data=relative_abundances,index=species_names).to_csv(output.TSV, sep='\t', header=False)
        cleaned_tax = list(zip(genus_names, species_names, relative_abundances))
        pd.DataFrame(cleaned_tax).to_csv(output.TAX_TSV, sep='\t', index=False, header=False)

rule sylph_output:
    input:
        cleaned_abundance_tax = rules.sylph_clean.output.TAX_TSV,
        classification = rules.sylph_taxprof.output.TSV

    output:
        TSV_sylph = ROOT / 'output' / '{class_level}' / 'output_sylph.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_sylph.log'
    run:
        import pandas as pd

        with open(input.classification, 'r') as f:
            for line in f:
                if line.split('\t')[0] == "UNKNOWN":
                    unknown_portion = 0
                    continue
                else:
                    unknown_portion = 0

        sylph_input = pd.read_table(input.cleaned_abundance_tax,
                                        names=['genus', 'species', 'abundance'])

        organism_count = pd.Series(sylph_input['abundance'].values / 100, index=sylph_input[params.level])
        organism_count['no_hit'] = unknown_portion / 100

        # Sum for the same organisms. This is important with genera as there can be multiple species with the same genus.
        organism_count = organism_count.groupby(organism_count.index).sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            braken2_input=organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV_sylph, sep='\t', index=True, header=False)
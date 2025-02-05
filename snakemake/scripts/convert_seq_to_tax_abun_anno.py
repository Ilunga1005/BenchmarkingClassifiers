import yaml
from pathlib import Path
import requests

s = requests.Session()


def fetch_reference_dataset_report(taxon, reference_only='true'):
    base_url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{taxon}/dataset_report?'
    params = {f'filters.reference_only': {reference_only},
              'filters.exclude_paired_reports': 'false',
              'filters.exclude_atypical': 'true'}
    try:
        response = s.get(base_url, params=params)
        # print(response.url)
        response.raise_for_status()
    except Exception as err:
        print(f'Other error occurred: {err}')
    else:
        return response.json()


def fetch_reference_sequence_report(accession):
    base_url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/sequence_reports?'
    params = {}
    try:
        response = s.get(base_url, params=params)
        response.raise_for_status()
    except Exception as err:
        print(f'Other error occurred: {err}')
    else:
        return response.json()


def fetch_fungus_or_nah(taxids):
    base_url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/'
    taxids_str = ','.join(map(str, taxids)) ## map()函数将taxids列表中的每个元素转换为字符串; ','.join()函数将字符串用，进行连结
    url = f'{base_url}{taxids_str}'
    params = {}
    try:
        response = s.get(url, params=params)
        response.raise_for_status() # 检查响应状态码，如果状态码不是200，则抛出HTTPError异常
    except Exception as err: # 异常处理
        print(f'Other error occurred: {err}')
    else: # 如果成功相应，则进入else模块
        r = response.json()

    fungi = {}
    for searched_item in r['taxonomy_nodes']:
        taxid = searched_item['taxonomy']['tax_id']
        if 4751 in searched_item['taxonomy']['lineage']:
            fungi[taxid] = True
        else:
            fungi[taxid] = False

    return fungi


def calculate_ti_from_si(ground_truth, lengths, fungi):
    """
    DOI: 10.1038/s41592-021-01141-3
    \frac{S_i}{T_i} = \frac{\sum_{i=0}^n \frac{R_i}{L_i}}{R_i} * L_i
    In the formula, R_i can be replaced by S_i because the sum of R_i is 1.
    This function only holds for species with a ploidy of 1.

    """
    sum_R = sum(ground_truth.values())
    sum_R_div_L = 0
    for taxid, R in ground_truth.items():
        P = 1
        if fungi[taxid]:
            P = 2
        sum_R_div_L += R / (P * lengths[taxid])
    # sum_R_div_L = sum([R / lengths[taxid] for taxid, R in ground_truth.items()])

    ti = {}
    for taxid, R in ground_truth.items():
        P = 1
        if fungi[taxid]:
            P = 2
        if R == 0:
            ti[taxid] = 0
        else:
            ti[taxid] = round((sum_R / sum_R_div_L) * (R / (P * lengths[taxid])), 6)

    return ti


def extract_accession_and_length(chosen_genome):
    accession = chosen_genome['accession']

    if (assembly_type := chosen_genome['assembly_info']['assembly_type']) == 'diploid': ## 检查组装类型是否为diploid；
        print(f"Found {assembly_type}")
        exit()
    lengths = 0
    r = fetch_reference_sequence_report(accession)
    sequences = {}
    for sequence in r['reports']:
        if sequence['assigned_molecule_location_type'] == 'Plasmid':
            pass
        elif sequence['assigned_molecule_location_type'] == 'Chromosome':
            sequences[sequence['refseq_accession']] = sequence['length']
            lengths += sequence['length']

    return sequences, lengths  # seqences是一个字典{accession: length};每个染色体的长度；lengths是整个基因组的长度


with open(snakemake.input[0], 'r') as f:
    theoretical_abundance = {'species': yaml.load(f, Loader=yaml.FullLoader)}

fungi = fetch_fungus_or_nah(list(theoretical_abundance['species'].keys())) ## 返回一个字典，{taxid: True/False}

reference_information, accessions, sequences, lengths, final_info = {}, {}, {}, {}, {}
for taxid, seq_abund in theoretical_abundance['species'].items():
    r = fetch_reference_dataset_report(taxid)
    if r:
        r = r['reports']
        if len(r) > 1:
            print(f'Multiple genomes found for {taxid}')
            print(f'{[i["accession"] for i in r]}')
        chosen_genome = r[0]
        accessions[taxid] = chosen_genome['accession']
        reference_information[taxid] = chosen_genome
        sequences[taxid], lengths[taxid] = extract_accession_and_length(chosen_genome)

    else:
        print(f'Nothing found for {taxid}, searching for non-reference')
        r = fetch_reference_dataset_report(taxid, reference_only='false')
        if not r.get('reports', None):
            print(f'Nothing found for {taxid}')
            break
        r = r['reports']
        if len(r) > 1:
            print(f'Multiple non-reference genomes found for {taxid}')
            level_to_search = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
            found_entries = []
            for level in level_to_search:
                matching_entries = [genome for genome in r
                                    if genome['assembly_info']['assembly_level'] == level]
                if matching_entries:
                    found_entries.extend(matching_entries)
                    break  # Exit the loop if entries are found
            r = found_entries
        chosen_genome = r[0]
        accessions[taxid] = chosen_genome['accession'] ## 染色体编号
        reference_information[taxid] = chosen_genome ## 选择的染色体信息
        sequences[taxid], lengths[taxid] = extract_accession_and_length(chosen_genome)

final_info = {taxid: {accessions[taxid]: {'sequences': sequences[taxid], 'final_length': lengths[taxid]}} for taxid in
              theoretical_abundance['species']}

tax_abun = calculate_ti_from_si(theoretical_abundance['species'], lengths, fungi) ## 返回理论丰度（ti）

with open(snakemake.output[0], 'w') as outfile:
    yaml.dump(tax_abun, outfile, default_flow_style=False, sort_keys=False)

with open(snakemake.log[0], 'w') as outfile:
    yaml.dump(final_info, outfile, default_flow_style=False, sort_keys=False)

from pyfaidx import Fasta

def get_genomic_seq(genome, coord_dict):
    """
    :param coord_dict: dict
    :return: dict {"seq":str DNA sequence}
    """
    coord_dict['start'] = int(coord_dict['start'])
    coord_dict['end'] = int(coord_dict['end'])
    for key in ['bp_upstream', 'bp_downstream']:
        if key in coord_dict.keys():
            coord_dict[key] = int(coord_dict[key])
        else:
            coord_dict[key] = 0

    # convert chromosome to match correct format, should work for 3 common cases: 1,Chr1,chr1
    coord_dict['chr'] = 'chr' + coord_dict['chr'][-1:]
    if coord_dict['orientation'] == '+':
        seq = genome[coord_dict['chr']][
              coord_dict['start'] - coord_dict['bp_upstream']: coord_dict['end'] + coord_dict['bp_downstream']]
    elif coord_dict['orientation'] == '-':
        seq = -genome[chr][
               coord_dict['start'] - coord_dict['bp_downstream']: coord_dict['end'] + coord_dict['bp_upstream']]
    else:
        raise AttributeError('Gene Orientation must be + or -!')
    return seq.seq


def get_gene_coords_noS3(event):
    f = open('TAIR10_genes.tab', 'r')
    storage_dict = {}
    # build a dictionary of the start, end, orientation
    for line in f.readlines():
        s_line = line.strip('\n').split('\t')
        if s_line[3] == event['gene_id'].upper():
            storage_dict[s_line[3]] = {'chr': s_line[0], 'start': int(s_line[1]), 'end': int(s_line[2]),
                                       'orientation': s_line[4]}
    # query dictionary for a specific gene
    try:
        result = storage_dict.get(event['gene_id'].upper())
    except KeyError:
        result = 'No match found'

    return result
    # raise Exception('Something went wrong')


def lambda_handler(event, context):
    gene_info = get_gene_coords_noS3(event)
    gene_info['name'] = event['gene_id']
    genome = Fasta('AT9.fa')
    output = get_genomic_seq(genome, gene_info)
    gene_info['sequence'] = output
    return gene_info
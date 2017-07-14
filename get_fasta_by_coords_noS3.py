from pyfaidx import Fasta
import json
import boto3

"""
{
    "q_info": {
    "chr": "$input.params('chr')",
    "start": "$input.params('start')",
    "end": "$input.params('end')",
    "orientation": "$input.params('orientation')"
    "bp_upstream": "$input.params('bp_upstream')"
    "bp_downstream": "$input.params('bp_downstream')"}
}
sample for testing
{
    "q_info": {
    "chr": "1",
    "start": "50",
    "end": "1000,
    "orientation": "+"
    "bp_upstream": "0"
    "bp_downstream": "0"}
}

"""

def get_genomic_seq(genome, coord_dict):
    """
    :param coord_dict: dict
    :return: dict {"seq":str DNA sequence}
    """
    coord_dict['start'] = int(coord_dict['start'])
    coord_dict['end'] = int(coord_dict['end'])
    coord_dict['bp_upstream'] = int(coord_dict['bp_upstream'])
    coord_dict['bp_downstream'] = int(coord_dict['bp_downstream'])
    # convert chromosome to match correct format, should work for 3 common cases: 1,Chr1,chr1
    coord_dict['chr'] = 'chr' + coord_dict['chr'][-1:]
    if coord_dict['orientation'] == '+':
        seq = genome[coord_dict['chr']][coord_dict['start'] - coord_dict['bp_upstream']: coord_dict['end'] + coord_dict['bp_downstream']]
    elif coord_dict['orientation'] == '-':
        seq = -genome[chr][coord_dict['start'] - coord_dict['bp_downstream']: coord_dict['end'] + coord_dict['bp_upstream']]
    else:
        raise AttributeError('Gene Orientation must be + or -!')
    return seq.seq


def lambda_handler(event, context):
    genome = Fasta('AT9.fa')
    output = get_genomic_seq(genome, event['q_info'])
    event['q_info']['sequence'] = output
    return event['q_info']

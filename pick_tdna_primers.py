import click
from pyfaidx import Fasta


def show_help_if_missing(ctx, param, value):
    if not value:
        click.secho("The '{0}' parameter is mandatory, please check the usage below:\n".format(param.name),
                    fg='yellow', bold=True)
        click.secho(ctx.get_help())
        ctx.exit()
    return value


@click.command()
@click.option('-i', '--input_stock_num', help='Stock number to generate primers for', callback=show_help_if_missing)
@click.option('-N', '--maxN', default=300, help='maxN or maximum distance from flanking sequence to T-DNA insertion')
@click.option('--ext5', default=300, help='ext5, basepairs upstream of T-DNA to exclude')
@click.option('--ext3', default=300, help='ext3, basepairs downstream of T-DNA to exclude')
@click.option('-z', '--p_zone', default=200, help='p_zone, size of region for picking each primer')
@click.option('-s', '--show-seq', default=False, help='Show the sequence around the insertion site', is_flag=True)
def run_tdna_primers(input_stock_num, maxn, ext5, ext3, p_zone, show_seq):
    # TODO: add option to show primer3 output details and save to file

    sequence_length_defs = {
        'maxN': maxn,
        'ext5': ext5,
        'ext3': ext3,
        'p_zone': p_zone * 2  # size of left and right regions together, so twice the Signal.salk.edu value
    }

    data_filenames = ['data/T-DNA.SALK', 'data/T-DNA.SAIL', 'data/T-DNA.GABI']
    db = init_db(data_filenames)
    genome = Fasta('data/AT9.fa')

    result = {}
    for poly_name, poly_info in db[input_stock_num].items():
        print('Making primers for ' + poly_name + '\n')
        sequence = get_seq(poly_info, sequence_length_defs, genome)
        if show_seq is True:
            click.echo('Genomic sequence around insert:\n')
            click.echo(sequence, color='red')
            click.echo('\n')
        result['poly_name'] = poly_name
        p3s, p3p = load_conf()
        primer_results = make_primers(sequence, sequence_length_defs, p3s, p3p)
        print('LP: ' + primer_results['PRIMER_LEFT_0_SEQUENCE'] + ' Tm: ' + str(
            round(primer_results['PRIMER_LEFT_0_TM'], 1)))
        print('RP: ' + primer_results['PRIMER_RIGHT_0_SEQUENCE'] + ' Tm: ' + str(
            round(primer_results['PRIMER_RIGHT_0_TM'], 1)))
        print('\n')


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def init_db(filenames):
    assert type(filenames) == list, "Filenames must be in the form of a list"
    records = {}
    for fname in filenames:
        f = open(fname, 'r')
        for line in f:
            s_line = line.strip('\n').split('\t')
            stock_name = s_line[0].split('.')[0]
            poly_name = s_line[0]
            poly_chr = s_line[1].split(':')[0]
            poly_chr = poly_chr[3:]
            poly_locs = s_line[3].split(',')[0]
            poly_start = poly_locs.split('/')[1].split('-')[0]
            poly_end = poly_locs.split('/')[1].split('-')[1]
            orientation = s_line[3].split('/')[0]
            # load into dict of dicts object
            if stock_name in records and poly_name in records[stock_name]:
                records[stock_name][poly_name + '.1'] = {'chr': poly_chr, 'orientation': orientation,
                                                         'start': poly_start,
                                                         'end': poly_end}
            elif stock_name not in records:
                records[stock_name] = {
                    poly_name: {'chr': poly_chr, 'orientation': orientation, 'start': poly_start, 'end': poly_end}}
            else:
                records[stock_name][poly_name] = {'chr': poly_chr, 'orientation': orientation, 'start': poly_start,
                                                  'end': poly_end}
    return records


def load_conf(conf_file='primer3.conf'):
    """
    loads all the settings for primer3 from the default file
    :param conf_file: string path to conf file
    :return: p3s: dict, p3p: dict
    """
    f = open(conf_file, 'r')
    p3s = {}
    p3p = {}
    for line in f.readlines():
        line = line.strip()
        value = line.split(',')[1]
        # divide into two separate dictionaries
        if line.startswith('SEQUENCE'):
            p3s[line.split(',')[0]] = value
        elif line.startswith('PRIMER'):
            # perform type conversions if needed
            # TODO: add handling for list arguments like PRIMER_PRODUCT_SIZE_RANGE
            if is_int(value):
                value = int(value)
            elif is_float(value):
                value = float(value)
            else:
                pass
            p3p[line.split(',')[0]] = value
        else:
            pass
    return p3s, p3p


def get_seq(poly_entry, sequence_length_defs, genome):
    """
    :param poly_entry: dict item
    'SALK_001127.52.15.x': {'chr': '1', 'orientation': 'W', 'start': '32', 'end': '303'}
    :param sequence_length_defs: dict
    :param genome: pyfaidx.Fasta
    :return seq: DNA sequence as string
    """
    # calculate the upstream and downstream distances based on the sequence input params
    bp_upstream = sequence_length_defs['maxN'] + sequence_length_defs['ext5'] + int(sequence_length_defs['p_zone'] / 2)
    bp_downstream = sequence_length_defs['ext3'] + int(sequence_length_defs['p_zone'] / 2)
    total_bp = bp_upstream + bp_downstream
    chrom = 'chr' + poly_entry['chr']

    # logic to account for orientation of T-DNA insert
    if poly_entry['orientation'] == 'W':
        new_start = max(1, int(poly_entry['start']) - bp_upstream)
        new_end = max(min(int(poly_entry['start']) + bp_downstream + 1, len(genome[chrom])), total_bp + new_start)
        seq = genome[chrom][new_start:new_end]
    elif poly_entry['orientation'] == 'C':
        print("reverse")
        new_start = max(int(poly_entry['end']) - bp_downstream, 1)
        new_end = max(min(int(poly_entry['end']) + bp_upstream, len(genome[chrom])), new_start + total_bp)
        seq = -genome[chrom][new_start:new_end]
    else:
        raise ValueError("Orientation must be 'W' or 'C'!")

    return str(seq)


def make_primers(sequence, sequence_length_defs, primer3_seq_args, primer3_primer_args):
    from primer3 import designPrimers
    from functools import reduce
    '''default iSECT values:
    size: optimal=21 min=18 max=28
    Tm: opt=61 Min=53 Max=71
    %GC: min=20 max=80
    clamp=1
    maxN=300 ext5=300 ext3=300
    primer_zone=200
    BPos=110 (distance from LB primer to insertion site)
    '''
    total = reduce(lambda x, y: x + y, [v for v in sequence_length_defs.values()])
    primer3_seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [0, int(sequence_length_defs['p_zone'] / 2),
                                                               total - int((sequence_length_defs['p_zone'] / 2)),
                                                               int((sequence_length_defs['p_zone'] / 2))]
    primer3_primer_args['PRIMER_PRODUCT_SIZE_RANGE'] = [total - sequence_length_defs['p_zone'], total]
    primer3_seq_args['SEQUENCE_TEMPLATE'] = sequence
    a = designPrimers(primer3_seq_args, primer3_primer_args)
    return a


if __name__ == '__main__':
    run_tdna_primers()

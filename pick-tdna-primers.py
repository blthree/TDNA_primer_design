import click
from pyfaidx import Fasta
from primer3.bindings import designPrimers


from tdna.sequence import get_seq
from tdna.IO import init_db, load_conf
from tdna.pick_primers import make_primers

def show_help_if_missing(ctx, param, value):
    if not value:
        click.secho("The '{0}' parameter is mandatory, please check the usage below:\n".format(param.name),
                    fg='yellow', bold=True)
        click.secho(ctx.get_help())
        ctx.exit()
    return value

@click.command()
@click.option('-i', '--input_stock_num', help='Stock number to generate primers for', callback=show_help_if_missing)
@click.option('-N','--maxN', default=300, help='maxN or maximum distance from flanking sequence to T-DNA insertion')
@click.option('--ext5', default=300, help='ext5, basepairs upstream of T-DNA to exclude')
@click.option('--ext3', default=300, help='ext3, basepairs downstream of T-DNA to exclude')
@click.option('-z', '--p_zone', default=200, help='p_zone, size of region for picking each primer')
def hello(input_stock_num, maxn, ext5, ext3, p_zone):
    sequence_length_defs = {
        'maxN': maxn,
        'ext5': ext5,
        'ext3': ext3,
        'p_zone': p_zone*2  # size of left and right regions together, so twice the Signal.salk.edu value
    }

    data_filenames = ['data/T-DNA.SALK', 'data/T-DNA.SAIL', 'data/T-DNA.GABI']
    db = init_db(data_filenames)
    genome = Fasta('data/AT9.fa')

    result = {}
    for poly_name, poly_info in db[input_stock_num].items():
        print('Making primers for ' + poly_name)
        sequence = get_seq(poly_info, sequence_length_defs, genome)
        result['poly_name'] = poly_name
        p3s, p3p = load_conf()
        primer_results = make_primers(sequence, sequence_length_defs, p3s, p3p)
        print('LP: ' + primer_results['PRIMER_LEFT_0_SEQUENCE'] + 'Tm: ' + str(round(primer_results['PRIMER_LEFT_0_TM'], 1)))
        #print('\n')
        print('RP: ' + primer_results['PRIMER_RIGHT_0_SEQUENCE'] + 'Tm: ' + str(round(primer_results['PRIMER_RIGHT_0_TM'], 1)))
        print('\n')
if __name__ == '__main__':
    hello()
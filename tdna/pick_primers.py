
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
    primer3_seq_args['p_zone'] = [0, int(sequence_length_defs['p_zone'] / 2),
                                  total - int((sequence_length_defs['p_zone'] / 2)), total]

    primer3_seq_args['SEQUENCE_TEMPLATE'] = sequence
    a = designPrimers(primer3_seq_args, primer3_primer_args)
    return a
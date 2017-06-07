from pyfaidx import Fasta
from primer3.bindings import designPrimers
from functools import reduce
data_filenames = ['data/T-DNA.SALK', 'data/T-DNA.SAIL', 'data/T-DNA.GABI']
# http://signal.salk.edu/database/transcriptome/AT9.fa
# http://signal.salk.edu/database/transcriptome/T-DNA.SALK
# The +/+ strand is 'W' on SALK data and +/- is 'C'
genome = Fasta('data/AT9.fa')
chrom_sizes = {
    'chr1': 30427671,
    'chr2': 19698289,
    'chr3': 23459830,
    'chr4': 18585056,
    'chr5': 26975502,
    'chrM':	366924,
    'chrC':	154478
}

primer3_seq_args = {
    'SEQUENCE_ID': 'T-DNA',
    'SEQUENCE_TEMPLATE': '',
    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0, 200, 1100, 200]
}
primer3_primer_args = {
    'PRIMER_TASK': 'generic',
    'PRIMER_LIBERAL_BASE': 1,
    'PRIMER_NUM_RETURN': 1,
    'PRIMER_EXPLAIN_FLAG': 0,
    'PRIMER_OPT_SIZE': 21,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 28,
    'PRIMER_OPT_TM': 56.5,
    'PRIMER_MIN_TM': 51.0,
    'PRIMER_MAX_TM': 61.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PRODUCT_SIZE_RANGE': [900, 1300],
    'PRIMER_GC_CLAMP': 1,
    'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 46
}
sequence_length_defs = {
    'maxN': 300,
    'ext5': 300,
    'ext3': 300,
    'p_zone': 400  # size of left and right regions together, so twice the Signal.salk.edu value
}
# calculate regions to pick primers from
total = reduce(lambda x, y: x + y, [v for v in sequence_length_defs.values()])
primer3_seq_args['p_zone'] = [0, int(sequence_length_defs['p_zone']/2), total-int((sequence_length_defs['p_zone']/2)), total]

def load_data(filenames):
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
                records[stock_name][poly_name + '.1'] = {'chr': poly_chr, 'orientation': orientation, 'start': poly_start,
                                                  'end': poly_end}
            elif stock_name not in records:
                records[stock_name] = {
                    poly_name: {'chr': poly_chr, 'orientation': orientation, 'start': poly_start, 'end': poly_end}}
            else:
                records[stock_name][poly_name] = {'chr': poly_chr, 'orientation': orientation, 'start': poly_start,
                                                  'end': poly_end}
    return records


def make_primers(sequence):
    '''default iSECT values:
    size: optimal=21 min=18 max=28
    Tm: opt=61 Min=53 Max=71
    %GC: min=20 max=80
    clamp=1
    maxN=300 ext5=300 ext3=300
    primer_zone=200
    BPos=110 (distance from LB primer to insertion site)
    '''
    primer3_seq_args['SEQUENCE_TEMPLATE'] = sequence
    a = designPrimers(primer3_seq_args, primer3_primer_args)
    return a

def get_seq(poly_entry):
    """
    :param poly_entry: dict item
    'SALK_001127.52.15.x': {'chr': '1', 'orientation': 'W', 'start': '32', 'end': '303'}
    :return seq: pyfaidx object
    """
    # calculate the upstream and downstream distances based on the sequence input params
    bp_upstream = sequence_length_defs['maxN'] + sequence_length_defs['ext5'] + sequence_length_defs['p_zone']
    bp_downstream = sequence_length_defs['ext3'] + sequence_length_defs['p_zone']
    total_bp = bp_upstream + bp_downstream
    chrom = 'chr' + poly_entry['chr']
    # logic to account for orientation of T-DNA insert
    print('Getting sequence from ' + chrom) 
    if poly_entry['orientation'] == 'W':
        # need to handle overflow cases, should define chromsizes!
        new_start = max(1, int(poly_entry['start']) - bp_upstream)
        new_end = max(min(int(poly_entry['start']) + bp_downstream + 1, len(genome[chrom])), total_bp + new_start)
        seq = genome[chrom][new_start:new_end]
    elif poly_entry['orientation'] == 'C':
        new_start = max(1, int(poly_entry['end']) - bp_downstream)
        new_end = max(min(int(poly_entry['end']) + bp_upstream + 1, len(genome[chrom])), new_start + total_bp)
        seq = -genome[chrom][new_start:new_end]
    else:
        raise ValueError("Orientation must be 'W' or 'C'!")
    print(str(seq))
    return str(seq)

def batch_process(db, out_file, number_to_make=1000, batch_size=100):
    f = open(out_file, 'a')
    primer_write_buffer = []
    for i in range(0, number_to_make):
        stock_num = list(db.keys())[i]
        for poly_name, poly_info in db[stock_num].items():
            print('Making primers for ' + poly_name)
            result = make_primers(get_seq(poly_info))
            result['poly_name'] = poly_name
            primer_write_buffer.append(result)
        # write to file and clear buffer after every batch_size items
        if i % batch_size == 0:
            if i != 0:
                for primer_pair in primer_write_buffer:
                    f.write("\t".join([str(v) for v in primer_pair.values()])+'\n')
                primer_write_buffer = []
                print('Wrote ' + str(batch_size) + ' primers to file: ' + out_file)

db = load_data(data_filenames)
batch_process(db, 'primers.txt', number_to_make=len(db))


def check_bases(seq):
    # count occurence of each letter in the sequence
    # easily check for ambiguous IUPAC symbols
    # that primer3 doesn't like
    # make sure to set the PRIMER_LIBERAL_BASE argument of primer3 to 1
    letters = {}
    for l in seq:
        if l in letters.keys():
            letters[l] += 1
        else:
            letters[l] = 1
    print(letters)
    return letters

def get_seq(poly_entry, sequence_length_defs, genome):
    """
    :param poly_entry: dict item
    'SALK_001127.52.15.x': {'chr': '1', 'orientation': 'W', 'start': '32', 'end': '303'}
    :return seq: DNA sequence as string
    """
    # calculate the upstream and downstream distances based on the sequence input params
    bp_upstream = sequence_length_defs['maxN'] + sequence_length_defs['ext5'] + sequence_length_defs['p_zone']
    bp_downstream = sequence_length_defs['ext3'] + sequence_length_defs['p_zone']
    total_bp = bp_upstream + bp_downstream
    chrom = 'chr' + poly_entry['chr']

    # logic to account for orientation of T-DNA insert
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

    return str(seq)

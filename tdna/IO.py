
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
                records[stock_name][poly_name + '.1'] = {'chr': poly_chr, 'orientation': orientation, 'start': poly_start,
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

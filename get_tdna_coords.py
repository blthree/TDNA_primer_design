import json
#event ={'tdna':'SALK_133000'}
def lambda_handler(event, context):
    records = {}
    f = open('T-DNA.SALK', 'r')
    query = event['tdna']
    for line in f.readlines():
        s_line = line.strip('\n').split('\t')
        if s_line[0].split('.')[0] == query:
            stock_name = s_line[0].split('.')[0]
            poly_name = s_line[0]
            poly_chr = s_line[1].split(':')[0]
            poly_chr = poly_chr[3:]
            poly_locs = s_line[3].split(',')[0]
            poly_start = poly_locs.split('/')[1].split('-')[0]
            poly_end = poly_locs.split('/')[1].split('-')[1]
            orientation = s_line[3].split('/')[0]
        # load into dict of dicts object
            if stock_name in records and poly_name in records[stock_name]['T-DNA insertions']:
                records[stock_name]['T-DNA insertions'][poly_name + '.1'] = {'Coordinates':{'chr': poly_chr, 'orientation': orientation,
                                                         'start': poly_start,
                                                         'end': poly_end}}
            elif stock_name not in records:
                records[stock_name] = {'T-DNA insertions': {poly_name: {'Coordinates': {'chr': poly_chr, 'orientation': orientation, 'start': poly_start, 'end': poly_end}}}}
            else:
                records[stock_name]['T-DNA insertions'][poly_name] = {'Coordinates':{'chr': poly_chr, 'orientation': orientation, 'start': poly_start,
                                                  'end': poly_end}}
    output = {'Stock Number': records}
    return records
#print(lambda_handler(event, ''))
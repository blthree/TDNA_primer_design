def parse():
    records = {}
    f = open('T-DNA.SALK', 'r')
    f_out = open('T-DNA.SALK.parsed', 'w')
    for line in f.readlines():
        s_line = line.strip('\n').split('\t')
        stock_name = s_line[0].split('.')[0]
        poly_name = s_line[0]
        poly_chr = s_line[1].split(':')[0]
        poly_chr = poly_chr[3:]
        poly_locs = s_line[3].split(',')[0]
        poly_start = poly_locs.split('/')[1].split('-')[0]
        poly_end = poly_locs.split('/')[1].split('-')[1]
        orientation = s_line[3].split('/')[0]
        orientation = orientation.replace('W', '+')
        orientation = orientation.replace('C', '-')
        # reassemble line
        fields = [stock_name, poly_name, poly_chr, poly_start, poly_end, orientation]
        simple_line = '\t'.join(fields) + '\n'
        #print(simple_line)
        f_out.write(simple_line)
import pandas as pd
#f = open('T-DNA.SALK.parsed', 'r')
df = pd.read_csv('T-DNA.SALK.parsed', sep='\t', names=['stock', 'poly', 'start', 'end', 'orientation'])
df.set_index('stock', inplace=True)
df.sort_index(inplace=True)
df.to_csv('T-DNA.SALK.sorted', sep='\t', header=None)


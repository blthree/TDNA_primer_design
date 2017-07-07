import click
from pyfaidx import Fasta

def make_ref_file():
    # CODE USED TO TRIM DOWN BED FILE TO TAB FILE
    f = open('TAIR10_genes.bed')
    f_out = open('TAIR10_genes.tab', 'w')
    for line in f.readlines():

        line = line.strip('\n').split('\t')
        # stop once we reach the plastid genome
        if line[0] == 'chrC':
            break
        # skip secondary isoforms e.g. AT1G12345.1
        if line[3][-2:-1] == '.':
            continue
        cut_line = line[0:4] + [line[5]]
        formatted_line = '\t'.join(cut_line)
        f_out.write(formatted_line + '\n')
    f.close()
    f_out.close()

f = open('TAIR10_genes.tab', 'r')
storage_dict = {}
# build a dictionary of the start, end, orientation
for line in f.readlines():
    s_line = line.strip('\n').split('\t')
    storage_dict[s_line[3]] = {'chr': int(s_line[0][3:]), 'start': int(s_line[1]), 'end': int(s_line[2]), 'orientation': s_line[4]}

# query dictionary for a specific gene
print(storage_dict.get('AT1G02580'))
# remember to convert all strings to uppercase before searching
print(storage_dict.get('at1g02580'.upper()))

class gene(object):
    def __init__(self, AGI):
        self.AGI = AGI.upper()
        self.dict_entry = storage_dict.get(self.AGI)
        # unpack the dictionary into properties
        self.chr, self.start, self.end, self.orientation = self.dict_entry.values()

    def get_genomic_seq(self):
        pass

    def get_expanded_seq(self):
        pass

gene1 = gene('AT1G02580')

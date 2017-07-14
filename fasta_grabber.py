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
        line[0] = line[0][-1:]
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
    storage_dict[s_line[3]] = {'chr': s_line[0], 'start': int(s_line[1]), 'end': int(s_line[2]), 'orientation': s_line[4]}

# query dictionary for a specific gene
print(storage_dict.get('AT1G02730'))
# remember to convert all strings to uppercase before searching
print(storage_dict.get('at1g02580'.upper()))

class gene(object):
    def __init__(self, AGI):
        self.AGI = AGI.upper()
        self.dict_entry = storage_dict.get(self.AGI)
        # unpack the dictionary into properties
        self.chr, self.start, self.end, self.orientation = self.dict_entry.values()
        self.chr = 'chr' + self.chr[-1:]

    def get_genomic_seq(self, bp_upstream=0, bp_downstream=0):

        if self.orientation == '+':
            seq = genome[self.chr][self.start - bp_upstream: self.end + bp_downstream]
            print(self.chr)
            print(self.start - bp_upstream)
            print(self.end + bp_downstream)
        elif self.orientation == '-':
            seq = -genome[self.chr][self.start - bp_downstream: self.end + bp_upstream]
        else:
            raise AttributeError('Gene Orientation must be + or -!')
        return seq

def get_seq(chr, start, end, orientation, bp_upstream=0, bp_downstream=0):
    chr = 'chr' + chr[-1:]
    if orientation == '+':
        seq = genome[chr][start - bp_upstream: end + bp_downstream]
    elif orientation == '-':
        seq = -genome[chr][start - bp_downstream: end + bp_upstream]
    else:
        raise AttributeError('Gene Orientation must be + or -!')
    return seq

genome = Fasta('AT9.fa')

gene1 = gene('AT1G02580')
print(gene1.get_genomic_seq())
a = gene1.get_genomic_seq()
print(type(a))

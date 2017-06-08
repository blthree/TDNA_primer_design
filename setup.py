import os, sys
import urllib.request
from setuptools import setup, find_packages

def reporthook(blocknum, blocksize, totalsize):
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar,))


setup(
    name='pick-tdna-primers',
    version='0.1',
    py_modules=['pick-tdna-primers'],
    install_requires=[
        'Click',
        'primer3-py',
        'pyfaidx',
    ],
    entry_points='''
        [console_scripts]
        pick-tdna-primers=TDNA_primer_design.pick-tdna-primers:run_tdna_primers
    ''',
)

def makeFastaIndex():
    from pyfaidx import Fasta
    sys.stderr.write('Building fasta index...')
    Fasta('data/AT9.fa')
    sys.stderr.write(' Done\n')

signal_url = 'http://signal.salk.edu/database/transcriptome/'
data_files = ['T-DNA.SALK', 'T-DNA.SAIL', 'T-DNA.GABI', 'AT9.fa']
if not os.path.exists('data'):
    os.mkdir('data', mode=755)

for file in data_files:
    sys.stderr.write('Checking for ' + file + '...')
    if not os.path.exists('data/' + file):
        sys.stderr.write(' not found, downloading...' + '\n')
        url_path = signal_url + file
        print(url_path)
        urllib.request.urlretrieve(url_path, 'data/' + file, reporthook)
    else:
        sys.stderr.write(' OK' + '\n')

makeFastaIndex()

sys.stderr.write('Creating default Primer3 config file...')
primer3_seq_args = {
    'SEQUENCE_ID': 'T-DNA',
    'SEQUENCE_TEMPLATE': ''
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
    'PRIMER_GC_CLAMP': 1,
    'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 46
}
f = open('primer3.conf', 'w')
for args in (primer3_seq_args, primer3_primer_args):
    for k,v in args.items():
        f.write(','.join((k,str(v)))+'\n')
f.close()
sys.stderr.write(' Done\n')
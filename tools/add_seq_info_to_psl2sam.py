import sys
from string import maketrans
import pysam
sam=sys.stdin
transtab = maketrans('ACTG','TGAC')

name2seq = {}
fa=open("/home/public/blat_test/hg19_CYP2D6_clean.fa")
for line in fa:
    name=line.strip()[1:]
    seq = fa.next().strip()
    name2seq[name]=seq

fa.close()

for line in sam:
    vals = line.split("\t")
    assert vals[0] in name2seq
    if vals[1]=='0':
        vals[-3]=name2seq[vals[0]]
    else:
        assert vals[1]=='16'
        vals[-3]=name2seq[vals[0]][::-1].translate(transtab)
    vals[5] = 'S'.join(vals[5].split('H'))
    sys.stdout.write("\t".join(vals))
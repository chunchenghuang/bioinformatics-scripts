import pysam
import numpy as np
import pandas as pd
from Bio import SeqIO
import sys
#from collections import Counter

dp=sys.argv[1]
#dp='./MA085_10_DP.txt'
sample_name=dp.split('_')[0]
con=sys.argv[2]
#con='./MA085_consensus.fa'
bed=sys.argv[3]
#bed='./SMN.bed'

def dpbedstat(chr,start,end):
    pos=start+1
    depth_value = depth_gc.loc[str(chr)].loc[int(pos):int(end)]['depth']
    mean = float(np.mean(depth_value))
    sd = float(np.std(depth_value))
    return mean,sd

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def movingcount(interval, window_size):
    window = np.ones(int(window_size))
    return np.convolve(interval, window, 'same')

def gccount(yourseq):
    yourseq_upper=map(str.upper,yourseq)
    cgs={'C':1,'G':1,'S':1}
    return pd.Series(yourseq_upper).map(cgs).fillna(0).astype(int)

depth_file=pd.read_csv(dp,sep="\t",comment='#',header=None)
depth_header = ['chr', 'pos', 'depth']
depth_file.columns = depth_header
depth_file = depth_file.set_index(['chr','pos'])

with open(con) as fasta_file:
    bed_chr=list(dict.fromkeys(depth_file.index.get_level_values('chr')))
    allgc = np.array([])
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # generator
        if seq_record.id in bed_chr:
            chrseq=list(seq_record.seq)
            chrgc=movingcount(gccount(chrseq),70)
            chrregion=depth_file.loc[seq_record.id].index-1
            bedgc=chrgc[chrregion]
            allgc=np.append(allgc,[seq_record.id,chrregion+1,bedgc])
        else:
            continue

fasta_file.close()

allgc_df=pd.DataFrame(allgc.reshape(len(bed_chr),3),columns=['chr','pos','gc']).set_index('chr')

unnested_lst = []
for col in allgc_df:
    unnested_lst.append(allgc_df[col].apply(pd.Series).stack())

gcdf=pd.concat(unnested_lst, axis=1, keys=allgc_df.columns)

gcdf=gcdf.reset_index().drop(columns='level_1')
gcdf['pos']=gcdf['pos'].astype(int)
gcdf['gc']=gcdf['gc'].astype(int)
gcdf=gcdf.set_index(['chr','pos'])

depth_gc=pd.merge(depth_file,gcdf,on=['chr','pos'],how='inner')

#GC content normalisation
rdmd=np.median(depth_gc['depth'])
mingc=min(depth_gc['gc'])
maxgc=max(depth_gc['gc'])
gcdict={}
for count in range(mingc,maxgc+1):
    normfactor=rdmd/np.median(depth_gc['depth'][depth_gc['gc']==count]) #change to mean temp, was median
    gcdict[int(count)]=float(normfactor)

for number, value in gcdict.items():
    if value == np.inf :
        gcdict[number]=1

depth_gc['gc']=depth_gc['gc'].map(gcdict)
depth_gc['depth']=depth_gc['depth']*depth_gc['gc']
depth_gc=depth_gc.drop(columns='gc')

#moving average
depth_gc['depth']=movingaverage(depth_gc['depth'],70)

bed_file=pd.read_csv(bed,sep="\t",comment='#',header=None)
header = ['chr', 'start', 'end', 'name', 'score', 'strand']
bed_file.columns = header

#fetch column with numbers use df[df.columns[X]], or use df.iloc[:,columns]
#featch rows with numbers use df.iloc[rows,:] or df.iloc[rows,] or df.iloc[rows]

#for loop to run over rows in pandas 
#https://cmdlinetips.com/2018/12/how-to-loop-through-pandas-rows-or-how-to-iterate-over-pandas-rows/
final_ls=[]
for line in bed_file.itertuples():
    mean,sd=dpbedstat(line.chr,line.start,line.end)
    output=[line.chr,line.start,line.end,line.name,mean,sd]
    final_ls.append(output)

final_cols=['chr','start','end','name','mean','sd']
final_df = pd.DataFrame(final_ls, columns=final_cols)
#final_df.to_csv(sample_name+"_cnv.csv")
final_df.to_csv(sample_name+"_PGK1_cnv.csv")
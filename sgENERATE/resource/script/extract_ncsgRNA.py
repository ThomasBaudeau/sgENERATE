import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
import time
from Bio import pairwise2
from matplotlib_venn import venn2_unweighted
from matplotlib_venn import venn3_unweighted
from math import sqrt



def found_read_peri(peri,filename,output):
    fname=open(filename).readlines()[0].replace('\n','')
    newsam=open(output,'w')
    bamFP = pysam.AlignmentFile(peri, "rb")
    for read in bamFP:
        try:
            if not (read.tags[-1][1]=='gRNA' or read.tags[-2][1]=='gRNA'):
                    if read.tags[-1][0]=='XC':
                        newsam.write('@'+fname+'_'+read.qname+'_'+read.reference_name+'\n')
                        newsam.write(read.seq+'\n')
                    else:
                       pass
            else:
                pass                
        except:
            print('exception')
            newsam.write('@'+fname+'_'+read.qname+'_'+read.reference_name+'\n')
            newsam.write(read.seq+'\n')

    return

found_read_peri(snakemake.input['a'],snakemake.input['b'], snakemake.output[0])
# found_read_peri('../../../Expcompmappeur/BIOREAL/tuned/Periscope_mult/COV_periscope.bam')

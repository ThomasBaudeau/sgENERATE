import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
import time
from Bio import pairwise2
from matplotlib_venn import venn2_unweighted
from matplotlib_venn import venn3_unweighted
from math import sqrt

TODEL={'peri':{},'peri2':{}}     


def tosave(peri,peri2):
    
    # open output bam with the header we just got
    bamFP = pysam.AlignmentFile(peri, "rb")
    bamFP2 = pysam.AlignmentFile(peri2, "rb")
    bam_header1 = bamFP.header.copy().to_dict()
    outbamfileperi = pysam.AlignmentFile('output_periscope.bam', "wb", header=bam_header1)
    bam_header2 = bamFP2.header.copy().to_dict()
    outbamfileperim = pysam.AlignmentFile('output_periscopem.bam', "wb", header=bam_header2)
    for read in bamFP:
        try:
            if not (read.tags[-1][1]=='gRNA' or read.tags[-2][1]=='gRNA'):
                if read.qname in TODEL['peri']:
                    outbamfileperi.write(read)
            else:
                if read.qname in TODEL['peri']:
                    outbamfileperi.write(read)
        except:
            pass
    for read2 in bamFP2:
        try:
            if not (read2.tags[-1][1]=='gRNA' or read2.tags[-2][1]=='gRNA'):
                if read2.qname in TODEL['peri']:
                    outbamfileperim.write(read2)
            else:
                if read2.qname in TODEL['peri']:
                    outbamfileperim.write(read2)
        except:
            pass
    raise

def found_read_peri(peri,GPvein,key,lname):
    
    # aresult={}
    print(key)
    for k in key:
        GPvein[lname][k]=[]
        # aresult[k]=0
    # aresult['MN908947.3']=0
    # aresult['ORF7b']=0
    # aresult['ORF1ab']=0
    newsam=open('onlynoncanonical.fasta','w')
    GPvein[lname]['non_canonical']=[]
    bamFP = pysam.AlignmentFile(peri, "rb")
    for read in bamFP:
        try:
            if not (read.tags[-1][1]=='gRNA' or read.tags[-2][1]=='gRNA'):
                    # read.tags[-2][1] 
                    # read.tags[-1][1]
                    if read.tags[-1][0]=='XC':
                        GPvein[lname]['non_canonical'].append(read.qname)
                        # if read.qname not in ok.keys():
                        #     ok[read.qname]=0
                        #     aresult[read.reference_name]+=1
                        # if lname!='peri':
                        #     TODEL[lname][read.qname]=read
                        # else:
                        #     if read.qname not in TODEL['peri']:
                        #         TODEL[lname][read.qname]=read

                        #         print(pairwise2.align.localms('AACCAACTTTCGATCTCTTGTAGATCTGTTCT', read.seq, 2, -2, -10, -.1,score_only=True))
                        #         print(read.qname)
                        # newsam.write('@'+read.qname+'\n')
                        # newsam.write(read.seq+'\n')
                        #
                    else:
                        GPvein[lname][read.tags[-1][1]].append(read.qname)
                        if lname=='peri2':
                            TODEL[lname][read.qname]=read
                        else:
                            if read.qname not in TODEL['peri2']:
                                TODEL[lname][read.qname]=read
                        
        except:
            # GPvein[lname]['non_canonical'].append(read.qname)
            pass
            # if read.qname not in ok.keys():
            #     ok[read.qname]=0
            #     aresult[read.reference_name]+=1
            #     if lname!='peri':
            # name=read.qname
            # theseq=read.seq
            # newsam.write('@'+name+'\n'+theseq+'\n')
    # for read,x  in enumerate(GPvein['peri2']['E']):
    #     print(read,x.qname,x.is_read1)
    return


def main(file1,file2,inperi,filenovel,file1mult,inperi2,filenovelmult,output,nb,file3=None):

    sgcount={}
    sgname={}
    ttreads=int(open(nb,'r').readlines()[0])/4
    gpvein={'peri2':{},'GT':{},'peri':{}}
    print('loading file')
    f1=open(file1,'r').readlines()
    f1p=open(file1mult,'r').readlines()
    f2=open(file2,'r').readlines()
    f4=open(filenovel,'r').readlines()
    f4p=open(filenovelmult,'r').readlines()
    for lign in f1[1:]:
        l=lign.split(',')
        sgcount[l[3]]=[int(l[4]),0,0]#
        gpvein['GT'][l[1]]=[]
    gpvein['GT']['non_canonical']=[]
    for lign in f1p[1:]:
        l=lign.split(',')
        try:
            sgcount[l[1]][2]=int(l[5])+int(l[6])#+int(l[7])
        except:
            sgcount[l[1]]=[0,int(l[5])+int(l[6]),0]
    for li in f2:
        l=li.split('_')
        if l[0]=='>sgRNA':
            sgname[int(l[3])]=l[1]
    nvcount=0
    for y in f4[1:]:
        ly=y.split(',')
        nvcount+=int(ly[3])
    sgcount['non_canonical']=[nvcount,0,0]
    nvcount=0
    for y2 in f4p[1:]:
        ly2=y2.split(',')
        nvcount+=int(ly2[5])+int(ly2[6])
    sgcount['non_canonical'][2]=nvcount
    sgRna=[]
    if file3!=None:
        f3=open(file3,'r').readlines()
        for x in f3:
            res=x.split('\t')
            num=int(res[0].split('_')[1].replace('.fastq',''))
            sg=sgname[num]
            sgcount[sg][1]+=int(res[1])
            for i in res[2:]:
                gpvein['GT'][sg].append(i)
        real=False
    else:
        real=True
    ttsgRna=0
    ttsgRnaP=0
    result={'Ground_truth':[],'Periscope':[],'Periscope_multi':[]}
    found_read_peri(inperi2,gpvein,list(sgcount.keys()),'peri2')
    found_read_peri(inperi,gpvein,list(sgcount.keys()),'peri')
    
    tosave(inperi,inperi2)
    for key in sgcount.keys():
            sgRna.append(key)
            if len(sgcount[key])==3:
                result['Ground_truth'].append(sgcount[key][1])  
                ttsgRna+=sgcount[key][1]
                result['Periscope'].append(sgcount[key][0])
                ttsgRnaP+=sgcount[key][0]
                result['Periscope_multi'].append(sgcount[key][2])
            else:
                result['Periscope'].append(sgcount[key][0])
                ttsgRnaP+=sgcount[key][1]
                result['Ground_truth'].append(0)
                result['Periscope_multi'].append(0)
    if real:
        ttsgRna=ttsgRnaP 
    print('data loaded')
    plot_lsgrna(result,sgRna,gpvein,ttreads,ttsgRna,output,real)


def plot_lsgrna(result,sgRna,gpvein,ttreads,ttsgRna,output,real):
        colors={'Ground_truth':'lightgreen','Periscope_multi':'coral','Periscope':'cornflowerblue'}
        add=1
        if real:
            result.pop('Ground_truth')
            add=0
        x = np.arange(int(len(sgRna)))*(3+add)  # the label locations
        width = 1  # the width of the bars
        multiplier = 0 
        fig=plt.figure(figsize=(16,13))
        sfigs=fig.subfigures(2, 1,wspace=0.2,hspace=0.2)   
        (ax2,ax) = sfigs[0].subplots(1,2,gridspec_kw={'width_ratios': [1, 10]})
        # fig, (ax2,ax) = plt.subplots(1,2,figsize=(16, 6),gridspec_kw={'width_ratios': [1, 10]})
        # fig.subplots_adjust(bottom=0.1)
        for attribute, measurement in result.items():
            offset = width * multiplier
            rects = ax.bar(x + offset, measurement, width, label=attribute,color=colors[attribute])
            ax.bar_label(rects, padding=5,fontsize=6)
            multiplier += 1
    
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('nb sgRNA')
        ax.set_title('nb sgRNA vs expected nb of sgRNA')
        add=0.5
        if real:
            add=0
        
        ax.set_xticks(x+0.5+add , sgRna)
        
        ax.legend(loc='upper right')
        if real:
            otherlabels=['Total\nRead','Total\nsgRNA\n(Periscope)']
        else:
            otherlabels=['Total\nRead','Total\nsgRNA']
        width2=0.85
        for ct,a in enumerate([ttreads,ttsgRna]):
            rects2 = ax2.bar( otherlabels[ct],a , width2)
            ax2.bar_label(rects2, padding=4)
        ax2.set_title('Sample\'s proportion',pad=15)
        ax2.set_ylabel('nb reads\n')
        ax2.set_xticks(np.arange(2) ,otherlabels)
        ax2.set_yticks([i for i in range(0,int(ttreads*1.2),int((ttreads*1.2)/5))])

        a=int(sqrt(len(gpvein['peri2'].keys())))
        b=int(len(gpvein['peri2'].keys())/int(sqrt(len(gpvein['peri2'].keys()))))
        axs_bottom=sfigs[1].subplots(a,b)
        if real:
            for key,axo in zip(gpvein['peri2'].keys(),axs_bottom.flat):
                v=venn2_unweighted([set(gpvein['peri2'][key]),set(gpvein['peri'][key])],ax=axo,set_colors=('coral', 'cornflowerblue'), alpha = 0.7)  
                v.get_label_by_id('A').set_text('Periscope_multi')
                v.get_label_by_id('B').set_text('Periscope')
                for text in v.set_labels:
                    text.set_fontsize(8)
                axo.set_title(key)
        else:
            for key,axo in zip(gpvein['peri2'].keys(),axs_bottom.flat):
                v=venn3_unweighted([set(gpvein['peri2'][key]), set(gpvein['GT'][key]),set(gpvein['peri'][key])],ax=axo,set_colors=('coral','lightgreen' ,'cornflowerblue'), alpha = 0.7) 
                v.get_label_by_id('A').set_text('Periscope_multi')
                v.get_label_by_id('B').set_text('ground_truth') 
                v.get_label_by_id('C').set_text('periscope')
                for text in v.set_labels:
                    text.set_fontsize(8)
                axo.set_title(key)
        print('Saving plot')    
        fig.savefig(output,dpi=300,format='pdf')
        
#main('Periscope/COV_periscope_counts.csv',  'result/COV_multifastq.faa','Periscope/COV_periscope.bam', 'Periscope/COV_periscope_novel_counts.csv','Periscope_mult/COV_periscope_counts.csv','Periscope_mult/COV_periscope.bam','Periscope_mult/COV_periscope_novel_counts.csv','test.pdf','result/final_COV_proportion.txt') #, 


# try:            
#     main(snakemake.input['a'],snakemake.input['b'],snakemake.input['peri'],snakemake.input['d'],snakemake.input['a2'],snakemake.input['peri2'],snakemake.input['d2'],snakemake.output[0],snakemake.input['nbread'],file3=snakemake.input['c']) 
# except AttributeError:
#     main(snakemake.input['a'],snakemake.input['b'],snakemake.input['peri'],snakemake.input['d'],snakemake.input['a2'],snakemake.input['peri2'],snakemake.input['d2'],snakemake.output[0],snakemake.input['nbread'])

#main("Periscope/COV_periscope_counts.csv", "result/COV_multifastq.faa", "Periscope/COV_periscope.bam", "Periscope/COV_periscope_novel_counts.csv", "Periscope_mult/COV_periscope_counts.csv", "Periscope_mult/COV_periscope.bam", "Periscope_mult/COV_periscope_novel_counts.csv",'result/final_COV_SG_ER.pdf' ,"result/COV_nbread.txt")
a='../Result_longmini/'
b='../Result_shortReadbwa/'
main(a+"Periscope_mult/COV_periscope_counts.csv", "result/COV_multifastq.faa", a+"Periscope_mult/COV_periscope.bam", a+"Periscope_mult/COV_periscope_novel_counts.csv", b+"Periscope_mult/COV_periscope_counts.csv", b+"Periscope_mult/COV_periscope.bam", b+"Periscope_mult/COV_periscope_novel_counts.csv",'result/final_COV_SG_ER.pdf' ,"result/COV_nbread.txt")
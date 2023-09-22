import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
import time
from Bio import pairwise2
from matplotlib_venn import venn2_unweighted
from matplotlib_venn import venn3_unweighted
from math import sqrt


def search_reads(read,search):
    """
    given a pysam read object and a search string perform a localms alignment
    :param read: pysam read object
    :param search: DNA search string e.g. ATGTGCTTGATGC
    :return: dictionary containing the read_id, alignment score and the position of the read
    """
    align_score = pairwise2.align.localms(search, read.seq, 2, -2, -10, -.1,score_only=True)
    return align_score


def read_minimap(f1,resu,sgname,gpvein):
    bamFP = pysam.AlignmentFile(f1, "rb")
    ttread=0
    count=0
    unmapped=0
    missmappedFN=0
    missmappedFP=0
    for keys in resu.keys():
        gpvein['minimap'][keys]=[]
        gpvein['GT'][keys]=[]
    for read in bamFP:
        if read.is_unmapped:
            count+=1
            ttread+=1
            if int(read.qname.split('_')[0].replace('S','')) in sgname.keys():
                gpvein['GT'][sgname[int(read.qname.split('_')[0].replace('S',''))]].append(read.qname)
                unmapped+=1
        elif read.is_secondary or read.is_supplementary:
            continue
        else:
            ttread+=1
            if len(read.reference_name.split('_'))<3:
                if int(read.qname.split('_')[0].replace('S',''))==int(read.reference_name.split('_')[1]):
                    # search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'
                    #sc=search_reads(read,search)
                    #if sc>50:
                        #notsupose+=1
                    continue
                else:
                    if int(read.qname.split('_')[0].replace('S','')) in sgname.keys():
                        missmappedFN+=1
                        gpvein['GT'][sgname[int(read.qname.split('_')[0].replace('S',''))]].append(read.qname)
                        
                        
            else:
                if int(read.qname.split('_')[0].replace('S',''))==int(read.reference_name.split('_')[3]):
                    search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'
                    #sc=search_reads(read,search)
                    resu[read.reference_name.split('_')[1]][2]+=1
                    gpvein['minimap'][read.reference_name.split('_')[1]].append(read.qname)
                    gpvein['GT'][sgname[int(read.qname.split('_')[0].replace('S',''))]].append(read.qname)
                else:
                    resu[read.reference_name.split('_')[1]][2]+=1
                    gpvein['minimap'][read.reference_name.split('_')[1]].append(read.qname)
                    missmappedFP+=1
                    if int(read.qname.split('_')[0].replace('S','')) in sgname.keys():
                        gpvein['GT'][sgname[int(read.qname.split('_')[0].replace('S',''))]].append(read.qname)
                        
    return resu,unmapped,missmappedFP,missmappedFN,ttread

def read_minimap_real(f1,resu,gpvein):
    bamFP = pysam.AlignmentFile(f1, "rb")
    ttread=0
    count=0
    nb_missed=0
    verylow=0
    notsupose=0
    for k in resu.keys():
        gpvein['minimap'][k]=[]
    for read in bamFP:
        if read.is_unmapped:
            ttread+=1
        elif read.is_secondary or read.is_supplementary:
            pass
        else:
            ttread+=1
            if read.reference_name.split('_')[1] in resu.keys():
                resu[read.reference_name.split('_')[1]][2]+=1
                gpvein['minimap'][read.reference_name.split('_')[1]].append(read.qname)
    return resu,nb_missed,notsupose,verylow,ttread

# def main(file1,file2,inperi,filenovel,mini,output,file3=None):
#     sgcount={}
#     sgname={}
#     gpvein={'minimap':{},'GT':{},'peri':{}}
#     f1=open(file1,'r').readlines()
#     f2=open(file2,'r').readlines()
#     f3=open(file3,'r').readlines()
#     f4=open(filenovel,'r').readlines()
#     for lign in f1[1:]:
#         l=lign.split(',')
#         sgcount[l[1]]=[int(l[5])+int(l[6])+int(l[7]),0]
#     for li in f2:
#         l=li.split('_')
#         if l[0]=='>sgRNA':
#             sgname[int(l[3])]=l[1]
#     for x in f3:
#         res=x.split('\t')
#         num=int(res[0].split('_')[1].replace('.fastq',''))
#         sg=sgname[num]
#         sgcount[sg][1]+=int(res[1])
#     nvcount=0
#     for y in f4[1:]:
#         ly=y.split(',')
#         nvcount+=int(ly[5])+int(ly[6])
#     sgcount['novel']=[nvcount,0]
#     sgRna=[]
#     ttsgRna=0
#     ttsgRnaP=0
#     result={'Ground_truth':[],'Periscope':[],'Minimap2_multi_ref':[]}
#     found_read_peri(inperi,gpvein,list(sgcount.keys()))
#     try:
#         sgcount,nb_missed,notsupose,verylow,ttreads=read_minimap(mini,sgcount,sgname,gpvein)
#         print(nb_missed,notsupose,verylow)
#         for key in sgcount.keys():
#             sgRna.append(key)
#             if len(sgcount[key])==3:
#                 result['Ground_truth'].append(sgcount[key][1])  
#                 ttsgRna+=sgcount[key][1]
#                 result['Periscope'].append(sgcount[key][0])
#                 ttsgRnaP+=sgcount[key][0]
#                 result['Minimap2_multi_ref'].append(sgcount[key][2])
#             else:
#                 result['Periscope'].append(sgcount[key][0])
#                 ttsgRnaP+=sgcount[key][1]
#                 result['Ground_truth'].append(0)
#                 result['Minimap2_multi_ref'].append(0)

#         x = np.arange(len(sgRna))  # the label locations
#         width = 0.25  # the width of the bars
#         multiplier = 0
#         fig=plt.figure(figsize=(16,13))
#         sfigs=fig.subfigures(2, 1,wspace=0.2,hspace=0.2)
#         (ax2,ax) = sfigs[0].subplots(1,2,gridspec_kw={'width_ratios': [1, 10]})
#         # fig, (ax2,ax) = plt.subplots(1,2,figsize=(16, 6),gridspec_kw={'width_ratios': [1, 10]})
#         # fig.subplots_adjust(bottom=0.1)
#         for attribute, measurement in result.items():
#             offset = width * multiplier
#             rects = ax.bar(x + offset, measurement, width, label=attribute)
#             ax.bar_label(rects, padding=5)
#             multiplier += 1
    
#         # Add some text for labels, title and custom x-axis tick labels, etc.
#         ax.set_ylabel('nb sgRNA')
#         ax.set_title('nb sgRNA vs expected nb of sgRNA')
#         ax.set_xticks(x + width, sgRna)
#         ax.legend(loc='upper right')

#         otherlabels=['Total\nRead','Total\nsgRNA']
#         width2=0.65
#         for ct,a in enumerate([ttreads,ttsgRna]):
#             rects2 = ax2.bar( otherlabels[ct],a , width2)
#             ax2.bar_label(rects2, padding=5)
#         ax2.set_title('Sample\'s proportion')
#         ax2.set_ylabel('nb reads')
#         a=int(sqrt(len(gpvein['minimap'].keys())))
#         b=int(len(gpvein['minimap'].keys())/int(sqrt(len(gpvein['minimap'].keys()))))
#         axs_bottom=sfigs[1].subplots(a,b)
#         for key,axo in zip(gpvein['minimap'].keys(),axs_bottom.flat):
#             v=venn3_unweighted([set(gpvein['minimap'][key]), set(gpvein['GT'][key]),set(gpvein['peri'][key])],ax=axo )  
#             v.get_label_by_id('A').set_text('minimap2')
#             v.get_label_by_id('B').set_text('ground_truth') 
#             v.get_label_by_id('C').set_text('periscope')
#             for text in v.set_labels:
#                 text.set_fontsize(8)
#             axo.set_title(key)
            
#         fig.savefig(output,dpi=550.00,format='pdf')
#     except:
#         print('REAL DATA\n######################################################################################################################################################################')
#         print(sgcount)
#         sgcount,nb_missed,notsupose,verylow,ttreads=read_minimap_real(mini,sgcount,gpvein)
#         print(sgcount)
#         for key in sgcount.keys():
#             sgRna.append(key)
#             if len(sgcount[key])==3:
#                 result['Ground_truth'].append(sgcount[key][1])  
#                 ttsgRna+=sgcount[key][1]
#                 result['Periscope'].append(sgcount[key][0])
#                 ttsgRnaP+=sgcount[key][0]
#                 result['Minimap2_multi_ref'].append(sgcount[key][2])
#             else:
#                 result['Periscope'].append(sgcount[key][0])
#                 ttsgRnaP+=sgcount[key][0]
#                 result['Ground_truth'].append(0)
#                 result['Minimap2_multi_ref'].append(0)

#         x = np.arange(len(sgRna))  # the label locations
#         width = 0.25  # the width of the bars
#         multiplier = 0
#         fig=plt.figure(figsize=(19,10))
#         sfigs=fig.subfigures(2, 1,wspace=0.2,hspace=0.2)
#         (ax2,ax) = sfigs[0].subplots(1,2,gridspec_kw={'width_ratios': [1, 10]})
#         result.pop('Ground_truth')
#         for attribute, measurement in result.items():
#             offset = width * multiplier
#             rects = ax.bar(x + offset, measurement, width, label=attribute)
#             ax.bar_label(rects, padding=5)
#             multiplier += 1
        
#         # Add some text for labels, title and custom x-axis tick labels, etc.
#         ax.set_ylabel('nb sgRNA')
#         ax.set_title('nb sgRNA vs expected nb of sgRNA')
#         ax.set_xticks(x + width, sgRna)
#         ax.legend(loc='upper right')

#         otherlabels=['Total\nRead','Total\nsgRNA\n(Periscope)']
#         width2=0.65
#         for ct,a in enumerate([ttreads,ttsgRnaP]):

#             rects2 = ax2.bar( otherlabels[ct],a , width2)
#             ax2.bar_label(rects2, padding=5)
        
        

#         ax2.set_title('Sample\'s proportion',x=0.5,y=1.15)
#         ax2.set_ylabel('nb reads')
#         a=round(sqrt(len(gpvein['minimap'].keys())))
#         b=round(len(gpvein['minimap'].keys())/sqrt(len(gpvein['minimap'].keys())))
#         print(a,b)
#         axs_bottom=sfigs[1].subplots(a,b)
#         for key,axo in zip(gpvein['minimap'].keys(),axs_bottom.flat):
#             v=venn2_unweighted([set(gpvein['minimap'][key]),set(gpvein['peri'][key])],ax=axo )  
#             v.get_label_by_id('A').set_text('minimap2')
#             v.get_label_by_id('B').set_text('periscope')
#             for text in v.set_labels:
#                 text.set_fontsize(8)
#             axo.set_title(key)
#         fig.savefig(output,dpi=550.00,format='pdf')
        

def found_read_peri(peri,GPvein,key):
    for k in key:
        GPvein['peri'][k]=[]
    bamFP = pysam.AlignmentFile(peri, "rb")
    newsam=open('onlynovel.bam','w')
    for read in bamFP:
       if not (read.tags[-1][1]=='gRNA' or read.tags[-2][1]=='gRNA'):
            read.tags[-2][1] 
            read.tags[-1][1]
            if read.tags[-1][0]=='XC':
                name=read.qname
                theseq=read.seq
                newsam.write('@'+name+'\n'+theseq+'\n')
                pass
            else:
                GPvein['peri'][read.tags[-1][1]].append(read.qname)
                if read.tags[-2][1].split('_')[1]=='LLQ':
                    GPvein['peri'][read.tags[-1][1]].pop()
                    pass
    return GPvein


def main(file1,file2,inperi,filenovel,mini,output,file3=None):
    sgcount={}
    sgname={}
    gpvein={'minimap':{},'GT':{},'peri':{}}
    f1=open(file1,'r').readlines()
    f2=open(file2,'r').readlines()
    f4=open(filenovel,'r').readlines()
    for lign in f1[1:]:
        l=lign.split(',')
        sgcount[l[1]]=[int(l[5])+int(l[6]),0,0]#+int(l[7])
    for li in f2:
        l=li.split('_')
        if l[0]=='>sgRNA':
            sgname[int(l[3])]=l[1]
    nvcount=0
    for y in f4[1:]:
        ly=y.split(',')
        nvcount+=int(ly[5])+int(ly[6])
    sgcount['novel']=[nvcount,0]
    sgRna=[]
    if file3!=None:
        f3=open(file3,'r').readlines()
        for x in f3:
            res=x.split('\t')
            num=int(res[0].split('_')[1].replace('.fastq',''))
            sg=sgname[num]
            sgcount[sg][1]+=int(res[1])
        real=False
    else:
        real=True
    ttsgRna=0
    ttsgRnaP=0
    result={'Ground_truth':[],'Periscope':[],'Minimap2_multi_ref':[]}
    found_read_peri(inperi,gpvein,list(sgcount.keys()))
    if not real:
        sgcount,nb_missed,notsupose,verylow,ttreads=read_minimap(mini,sgcount,sgname,gpvein)
    else:
        sgcount,nb_missed,notsupose,verylow,ttreads=read_minimap_real(mini,sgcount,gpvein)
    for key in sgcount.keys():
            sgRna.append(key)
            if len(sgcount[key])==3:
                result['Ground_truth'].append(sgcount[key][1])  
                ttsgRna+=sgcount[key][1]
                result['Periscope'].append(sgcount[key][0])
                ttsgRnaP+=sgcount[key][0]
                result['Minimap2_multi_ref'].append(sgcount[key][2])
            else:
                result['Periscope'].append(sgcount[key][0])
                ttsgRnaP+=sgcount[key][1]
                result['Ground_truth'].append(0)
                result['Minimap2_multi_ref'].append(0)
    if real:
        ttsgRna=ttsgRnaP 
    plot_lsgrna(result,sgRna,gpvein,ttreads,ttsgRna,output,real)


def plot_lsgrna(result,sgRna,gpvein,ttreads,ttsgRna,output,real):
        colors={'Ground_truth':'lightgreen','Minimap2_multi_ref':'coral','Periscope':'cornflowerblue'}
        if real:
            result.pop('Ground_truth')
        x = np.arange(len(sgRna))  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0
        fig=plt.figure(figsize=(16,13))
        sfigs=fig.subfigures(2, 1,wspace=0.2,hspace=0.2)
        
        (ax2,ax) = sfigs[0].subplots(1,2,gridspec_kw={'width_ratios': [1, 10]})
        # fig, (ax2,ax) = plt.subplots(1,2,figsize=(16, 6),gridspec_kw={'width_ratios': [1, 10]})
        # fig.subplots_adjust(bottom=0.1)
        for attribute, measurement in result.items():
            offset = width * multiplier
            rects = ax.bar(x + offset, measurement, width, label=attribute,color=colors[attribute])
            ax.bar_label(rects, padding=5)
            multiplier += 1
    
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('nb sgRNA')
        ax.set_title('nb sgRNA vs expected nb of sgRNA')
        ax.set_xticks(x + width, sgRna)
        ax.legend(loc='upper right')
        if real:
            otherlabels=['Total\nRead','Total\nsgRNA\n(Periscope)']
        else:
            otherlabels=['Total\nRead','Total\nsgRNA']
        width2=0.65
        for ct,a in enumerate([ttreads,ttsgRna]):
            rects2 = ax2.bar( otherlabels[ct],a , width2)
            ax2.bar_label(rects2, padding=5)
        ax2.set_title('Sample\'s proportion')
        ax2.set_ylabel('nb reads')
        a=int(sqrt(len(gpvein['minimap'].keys())))
        b=int(len(gpvein['minimap'].keys())/int(sqrt(len(gpvein['minimap'].keys()))))
        axs_bottom=sfigs[1].subplots(a,b)
        if real:
            for key,axo in zip(gpvein['minimap'].keys(),axs_bottom.flat):
                v=venn2_unweighted([set(gpvein['minimap'][key]),set(gpvein['peri'][key])],ax=axo,set_colors=('coral', 'cornflowerblue'), alpha = 0.7)  
                v.get_label_by_id('A').set_text('minimap2')
                v.get_label_by_id('B').set_text('periscope')
                for text in v.set_labels:
                    text.set_fontsize(8)
                axo.set_title(key)
        else:
            for key,axo in zip(gpvein['minimap'].keys(),axs_bottom.flat):
                v=venn3_unweighted([set(gpvein['minimap'][key]), set(gpvein['GT'][key]),set(gpvein['peri'][key])],ax=axo,set_colors=('coral','lightgreen' ,'cornflowerblue'), alpha = 0.7) 
                v.get_label_by_id('A').set_text('minimap2')
                v.get_label_by_id('B').set_text('ground_truth') 
                v.get_label_by_id('C').set_text('periscope')
                for text in v.set_labels:
                    text.set_fontsize(8)
                axo.set_title(key)    
        fig.savefig(output,dpi=550.00,format='pdf')
#main('Periscope/COV_periscope_counts.csv',  'result/COV_multifastq.faa','Periscope/COV_periscope.bam', 'Periscope/COV_periscope_novel_counts.csv', 'result/minimap_COV.bam','test.pdf') #, 'result/final_COV_proportion.txt'


try:            
    main(snakemake.input['a'],snakemake.input['b'],snakemake.input['peri'],snakemake.input['d'],snakemake.input['mini'],snakemake.output[0],file3=snakemake.input['c']) 
except AttributeError:
    main(snakemake.input['a'],snakemake.input['b'],snakemake.input['peri'],snakemake.input['d'],snakemake.input['mini'],snakemake.output[0])



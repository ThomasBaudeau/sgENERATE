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

def found_read_peri(peri,GPvein,key,lname,LLQ=True):
    
    # aresult={}
    for k in key:
        GPvein[lname][k]=[]
        # aresult[k]=0
    # aresult['MN908947.3']=0
    # aresult['ORF7b']=0
    # aresult['ORF1ab']=0
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
                        #     if lname!='peri':
                        #         print(pairwise2.align.localms('AACCAACTTTCGATCTCTTGTAGATCTGTTCT', read.seq, 2, -2, -10, -.1,score_only=True))
                        #         print(read.qname)
                        # newsam.write('@'+read.qname+'\n')
                        # newsam.write(read.seq+'\n')
                        # print(read.reference_name)
                        TODEL[lname][read.qname]=read
                    else:
                        GPvein[lname][read.tags[-1][1]].append(read.qname)
                        if read.tags[-2][1].split('_')[1]=='LLQ' and not LLQ:
                            GPvein[lname][read.tags[-1][1]].pop()
                        
        except:
            GPvein[lname]['non_canonical'].append(read.qname)
            # if read.qname not in ok.keys():
            #     ok[read.qname]=0
            #     aresult[read.reference_name]+=1
            #     if lname!='peri':
            #         name=read.qname
            #         theseq=read.seq
            #         newsam.write('@'+name+'\n'+theseq+'\n')
    return


def extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb,file3=None,LLQ=True):
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
        if LLQ:
            sgcount[l[1]]=[int(l[5])+int(l[6])+int(l[7]),0,0]#
        else:
            sgcount[l[1]]=[int(l[5])+int(l[6]),0,0]
        gpvein['GT'][l[1]]=[]
    gpvein['GT']['non_canonical']=[]
    for lign in f1p[1:]:
        l=lign.split(',')
        sgcount[l[1]][2]=int(l[5])+int(l[6])#+int(l[7])
    for li in f2:
        l=li.split('_')
        if l[0]=='>sgRNA':
            sgname[int(l[3])]=l[1]
    nvcount=0
    for y in f4[1:]:
        ly=y.split(',')
        nvcount+=int(ly[5])+int(ly[6])
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
            res=x.split('\t')[:-2]
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
    return result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname

def main(file1,file2,inperi,filenovel,file1mult,inperi2,filenovelmult,output,nb,file3=None):

    if file3:
        result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname =extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb,file3)
    else:
        result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname =extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb)
    found_read_peri(inperi,gpvein,list(sgcount.keys()),'peri')
    found_read_peri(inperi2,gpvein,list(sgcount.keys()),'peri2') 
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
        
# main('Periscope/COV_periscope_counts.csv',  'result/COV_multifastq.faa','Periscope/COV_periscope.bam', 'Periscope/COV_periscope_novel_counts.csv','Periscope_mult/COV_periscope_counts.csv','Periscope_mult/COV_periscope.bam','Periscope_mult/COV_periscope_novel_counts.csv','test.pdf','result/final_COV_proportion.txt') #, 


# def calc_f1score(gpvein):
#     TP1=0
#     FP1=0
#     FN1=0
#     TP2=0
#     FP2=0
#     FN2=0
#     for key in gpvein['peri2'].keys():
#         if key=='N*':
#             continue
#         curTP1=len(set(gpvein['peri'][key])&set(gpvein['GT'][key]))
#         TP1+= len(set(gpvein['peri'][key])&set(gpvein['GT'][key]))
#         FP1+=len(set(gpvein['peri'][key]))-curTP1
#         FN1+=len(set(gpvein['GT'][key]))-curTP1
#         curTP2=len(set(gpvein['peri2'][key])&set(gpvein['GT'][key]))
#         TP2+= len(set(gpvein['peri2'][key])&set(gpvein['GT'][key]))
#         FP2+=len(set(gpvein['peri2'][key]))-curTP2
#         FN2+=len(set(gpvein['GT'][key]))-curTP2
#     F11=TP1/(TP1+(0.5*(FN1+FP1)))
#     F22=TP2/(TP2+(0.5*(FN2+FP2)))
#     print('peri:')
#     print(F11)
#     print('perim')
#     print(F22)
#     raise



# d='BIOSIMU2'
# a='../../../../Expcompmappeur/'+d+'/default/Periscope/'
# b='../../../../Expcompmappeur/'+d+'/default/Periscope_mult/'
# c='../../../../Expcompmappeur/'+d+'/default/result/'
# print(os.getcwd())
# main(a+"COV_periscope_counts.csv", c+"COV_multifastq.faa", a+"COV_periscope.bam", a+"COV_periscope_novel_counts.csv", b+"COV_periscope_counts.csv", b+"COV_periscope.bam", b+"COV_periscope_novel_counts.csv",c+'ERR4_withLLQ.pdf' ,c+"COV_nbread.txt",c+'final_COV_proportion.txt')

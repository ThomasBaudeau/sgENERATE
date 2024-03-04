#!/usr/bin/env python3
#from Bio import pairwise2
import pysam
from BCBio import GFF
from Bio import pairwise2
from pybedtools import *
import glob

import pprint as pp
from Bio import SeqIO
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor as ProcessPool
class PeriscopeRead(object):
    def __init__(self, read):
        self.read = read



def create_group_ref(ul):

    #from a list of read return a list of read by reference

    a={}
    for i in ul:
        if not i.reference_name:
            continue
        else:
            try:
                a[i.reference_name].append(i)
            except:
                a[i.reference_name]=[i]   
    return a                 
                   
def create_group_pos(ul):
    #from a list of read return a list of read by start position on reference
    a={}
    for i in ul:
        try:
            a[i.reference_start].append(i)
        except:
            a[i.reference_start]=[i]   
    return a 


def correct_position(gff,read):
    #correct position between multireference and normal reference
    st=int(gff[read.reference_name]['begin'])
    nd=int(gff[read.reference_name]['begin'])
    if read.reference_name=='ORF1ab':
        return read.reference_start,read.reference_end
    else:
        if read.reference_start<36:
            nd+=read.reference_end-36
            return st,nd
        else:
            st+=read.reference_start-36
            nd+=read.reference_end-36
            return st,nd

def load_gff(gff_file):
    # create dictionary of each real starting positon for each gene
    in_file = gff_file

    in_handle = open(in_file)
    ab= GFF.parse(in_handle)
    end=None
    begin=None
    gff={}
    for rec in ab:
        for feature in rec.features:
            if feature.type=='gene':
                end=feature.location.nofuzzy_end
                if begin !=None:
                    if begin>feature.location.nofuzzy_start:
                        begin=feature.location.nofuzzy_start
                    gff[str(feature.qualifiers['Name'][0])]={'begin':begin,'end':end}
                else:
                    gff[str(feature.qualifiers['Name'][0])]={'begin':0,'end':end}
                begin=feature.location.nofuzzy_end
    in_handle.close()
    return gff


def cut_leader():
    #from a fasta file return a two fasta file without the leader part for each read
    nfile2=open('newsgRNA.fasta','w')
    nfile=open('newsgRNApluslong.fasta','w')
    for l in open('all_ncsgRNA.fasta').readlines():
        if l=='\n' or l[0]=='@':
            a=l
        else:
            align= pairwise2.align.localms(l,'AACCAACTTTCGATCTCTTGTAGATCTGTTCT' ,2, -2, -10, -.1)
            reste=l[align[0].end:max([align[0].end+50,(len(l)-25)])]
            reste2=l[align[0].end:min([align[0].end+200,(len(l)-50)])]
            nfile.write(a.replace('\n','')+'\n')
            nfile.write(reste.replace('\n','')+'\n')
            nfile2.write(a.replace('\n','')+'\n')
            nfile2.write(reste2.replace('\n','')+'\n')


def is_truncated():
    #for each result of ORFfinder tell if the protein will be truncated or not 
    tsv=open('all_nc_sgrna_1.csv','w')
    tsv.write('sample;name;junction_position;gene;protein+\n')
    c1=0
    c2=0
    for f1 in glob.glob('alignProt/result/*'):
        a=True
        sample=f1.split('_')[0].replace('alignProt/result/','')
        name='_'.join(f1.split('_')[1:3])
        name=sample+'_'+name
        dec=f1.split('_')[-1].replace('.fasta','')
        ref=f1.split('_')[-2].replace('.fasta','')
        f2='alignProt/'+ref+'_fin.fasta'
        l1=[]
        l2=[]
        if int(dec)>40000:
            pass
        for seq_record in SeqIO.parse(f1, "fasta"):
            l1.append([int(seq_record.id.split(':')[1]),seq_record.seq[0:50]])
        for seq_record in SeqIO.parse(f2, "fasta"):
            l2.append([int(seq_record.id.split(':')[1]),seq_record.seq[0:50]])
            if seq_record.description.split(':')[-1]=='wait':
                waited=seq_record.seq[0:50]
                pos=seq_record.id.split(':')[1]
        l1=sorted(l1)
        for i in range(len(l1)):
            if l1[i][1]==waited:
                tsv.write(sample+';'+name+';'+dec+';'+ref+';'+'normal'+'\n')
                a=False
                if ref=='S':
                    print(':)')
        if a:
            tsv.write(sample+';'+name+';'+dec+';'+ref+';'+'truncated'+'\n')



def compare_4_alignment():
    #from 4 alignment extract a concensus reference position and return a sequence from this position
    inbamfile = pysam.AlignmentFile('all.bam', "rb")
    inbamfile1 = pysam.AlignmentFile('short.bam', "rb")
    inbamfile2 = pysam.AlignmentFile('long.bam', "rb")
    inbamfile3 = pysam.AlignmentFile('splice.bam', "rb")#pluslong
    dicref={}
    nbref={}
    av='/'
    for l in open('multi_ref.fa').readlines():
        if l[0]=='>':
            
            name=l.replace('\n','').replace('>','')
            nbref[av]=name
            av=name
        else:
            dicref[name]=l


    ref={}
    dicread={}
    counta=0
    countb=0
    countc=0
    countd=0
    gff='covid.gff'
    freste=open('some_nc_sgRNA.fasta','w')
    fresult=open('seq_ncsgRNA.fasta','w')
    gffa=load_gff(gff)
    tsv=open('all_nc_sgrna_3.csv','w')
    tsv.write('sample;name;junction_position;gene;protein+\n')
    for i in open('multi_ref.fa','r'):
        if i[0]=='>':
            name=i[1:].replace('\n','')
        else:
            ref[name]=i
    for read in inbamfile:
        if not (read.is_secondary or read.is_supplementary or read.is_unmapped):
            if read.reference_name!=read.qname.split('_')[-1]:
                pass
            dicread[read.qname]=[read]
    for read in inbamfile1:
        if not (read.is_secondary or read.is_supplementary ):
            try:
                dicread[read.qname].append(read)
            except:
                pass
    for read in inbamfile3:
        if not (read.is_secondary or read.is_supplementary ):
            try:
                dicread[read.qname].append(read)
            except:
                pass
    for read in inbamfile2:
        if not (read.is_secondary or read.is_supplementary ):
            try:
                dicread[read.qname].append(read)
                a=dicread[read.qname]
            except :
                continue
                # align0= pairwise2.align.localms(dicread[read.qname][0].seq,ref[dicread[read.qname][0].reference_name],2, -2, -10, -.1)
                # align1= pairwise2.align.localms(dicread[read.qname][1].seq,ref[dicread[read.qname][0].reference_name],2, -2, -10, -.1)
                # align2=pairwise2.align.localms(read.seq,ref[read.reference_name],2, -2, -10, -.1)
                # print()
                # print(align1[0].start,align2[0].start)
            try:
                st,nd=correct_position(gffa,a[0])
                
                ol=[a[0].reference_start,a[1].reference_start,a[2].reference_start,a[3].reference_start]
            except :
                continue
            pop=list(set(ol))
            
            if a[0].reference_name == a[1].reference_name ==a[2].reference_name==a[3].reference_name:
                if a[0].reference_name==read.qname.split('_')[-1] and len(pop)==1:
                    counta+=1
                    continue
                    
                    name='_'.join(read.qname.split('_')[:-1])
                    refn=read.qname.split('_')[-1]
                    fseq=dicref[read.qname.split('_')[-1]][pop[0]:]
                    change=False
                    if read.qname.split('_')[-1] !='ORF10':
                        change=True
                        seq=dicref[read.qname.split('_')[-1]][read.reference_start:gffa[nbref[read.qname.split('_')[-1]]]['begin']-st]
                        for i in range(len(seq)-3):
                            if seq[i]=='A' and seq[i+1]=='T' and seq[i+2]=='G':
                                change=False
                                break
                    if change:
                        fseq=dicref[nbref[read.qname.split('_')[-1]]]
                        refn=nbref[refn]
                    fresult.write('>'+name+'_'+refn+'_'+str(st)+'\n')
                    fresult.write(fseq)
                elif a[1].reference_start==a[2].reference_start==a[3].reference_start:
                    countb+=1
                    read=a[1]
                    continue
                    name='_'.join(read.qname.split('_')[:-1])
                    refn=read.qname.split('_')[-1]
                    fseq=dicref[read.qname.split('_')[-1]][a[3].reference_start:]
                    change=False
                    if read.qname.split('_')[-1] !='ORF10':
                        change=True
                        seq=dicref[read.qname.split('_')[-1]][a[3].reference_start:gffa[nbref[read.qname.split('_')[-1]]]['begin']-st]
                        change=True
                        for i in range(len(seq)-3):
                            if seq[i]=='A' and seq[i+1]=='T' and seq[i+2]=='G':
                                change=False
                                
                                break
                    if change:
                        fseq=dicref[nbref[read.qname.split('_')[-1]]]
                        refn=nbref[refn]
                    fresult.write('>'+name+'_'+refn+'_'+str(st)+'\n')
                    fresult.write(fseq)
                else:
                    countd+=1
                    #continue
                    a=dicread[read.qname]
                    sample=a[0].qname.split('_')[0]
                    name_read=a[0].qname
                    name=sample+'_'+name
                    st,nd=correct_position(gffa,a[0])
                    dec=str(st)
                    ref_name=a[0].reference_name
                    tsv.write(sample+';'+name_read+';'+dec+';'+ref_name+';'+'uncertain'+'\n')   
                                                     
            else:
                all_read=[a[0],a[1],a[3],a[2]]
                ul= [a[0].reference_name,a[1].reference_name,a[3].reference_name,a[2].reference_name]
                pup=list(set(ul))
                if a[0].reference_name == a[1].reference_name ==a[3].reference_name:
                    ol=[a[0].reference_start ,a[1].reference_start,a[3].reference_start]
                    pop=list(set(ol))
                    if len(pop)==1:
                        countb+=1
                        continue
                        read=a[0]
                        st,nd=correct_position(gffa,read)
                        name='_'.join(read.qname.split('_')[:-1])
                        refn=read.qname.split('_')[-1]
                        fseq=dicref[read.qname.split('_')[-1]][read.reference_start:]
                        change=False
                        if read.qname.split('_')[-1] !='ORF10':
                            seq=dicref[read.qname.split('_')[-1]][read.reference_start:gffa[nbref[read.qname.split('_')[-1]]]['begin']-st]
                            change=True
                            for i in range(len(seq)-3):
                                if seq[i]=='A' and seq[i+1]=='T' and seq[i+2]=='G':
                                    change=False
                                    break
                        if change:
                            fseq=dicref[nbref[read.qname.split('_')[-1]]]
                            refn=nbref[refn]
                        fresult.write('>'+name+'_'+refn+'_'+str(st)+'\n')
                        fresult.write(fseq)
                    else:
                        countd+=1
                        #continue
                        a=dicread[read.qname]
                        sample=a[0].qname.split('_')[0]
                        name_read=a[0].qname
                        name=sample+'_'+name
                        st,nd=correct_position(gffa,a[0])
                        dec=str(st)
                        ref_name=a[0].reference_name
                        tsv.write(sample+';'+name_read+';'+dec+';'+ref_name+';'+'uncertain'+'\n')
                elif len(pup)>=2:
                    gp_ref=create_group_ref(all_read)
                    ls_pos={}
                    for key in gp_ref.keys():
                        ls_pos[key]=[]
                        for read in gp_ref[key]:
                            ls_pos[key].append(read.reference_start)
                            ok=True
                    for key in ls_pos.keys():
                        if len(list(set(ls_pos[key])))!=1:
                            ok=False
                    if ok:
                        if len(gp_ref.keys())!=1:
                            countc+=1
                            #continue
                            a=dicread[read.qname]
                            sample=a[0].qname.split('_')[0]
                            name_read=a[0].qname
                            name=sample+'_'+name
                            st,nd=correct_position(gffa,a[0])
                            dec=str(st)
                            ref_name=a[0].reference_name
                            tsv.write(sample+';'+name_read+';'+dec+';'+ref_name+';'+'uncertain'+'\n')
                        else:
                            countb+=1
                            continue
                            read=gp_ref[list(gp_ref.keys())[0]][0]
                            st,nd=correct_position(gffa,read)
                            name='_'.join(read.qname.split('_')[:-1])
                            refn=read.qname.split('_')[-1]
                            fseq=dicref[read.qname.split('_')[-1]][read.reference_start:]
                            change=False
                            if read.qname.split('_')[-1] !='ORF10':
                                seq=dicref[read.qname.split('_')[-1]][read.reference_start:gffa[nbref[read.qname.split('_')[-1]]]['begin']-st]
                                change=True
                                for i in range(len(seq)-3):
                                    if seq[i]=='A' and seq[i+1]=='T' and seq[i+2]=='G':
                                        change=False
                                        
                                        break
                            if change:
                                fseq=dicref[nbref[read.qname.split('_')[-1]]]
                                refn=nbref[refn]
                            fresult.write('>'+name+'_'+refn+'_'+str(st)+'\n')
                            fresult.write(fseq)
                    else:
                        pass
                        countd+=1
                        #continue
                        a=dicread[read.qname]
                        sample=a[0].qname.split('_')[0]
                        name_read=a[0].qname
                        name=sample+'_'+name
                        st,nd=correct_position(gffa,a[0])
                        dec=str(st)
                        ref_name=a[0].reference_name
                        tsv.write(sample+';'+name_read+';'+dec+';'+ref_name+';'+'uncertain'+'\n')
                else:
                    
                    countd+=1
                    #continue
                    a=dicread[read.qname]
                    sample=a[0].qname.split('_')[0]
                    name_read=a[0].qname
                    name=sample+'_'+name
                    st,nd=correct_position(gffa,a[0])
                    dec=str(st)
                    ref_name=a[0].reference_name
                    tsv.write(sample+';'+name_read+';'+dec+';'+ref_name+';'+'uncertain'+'\n')
                                    
    print(counta,countb,countc,countd,counta+countb+countc+countd)         


                   
          
    
 




            
            

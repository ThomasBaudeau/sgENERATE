import pysam
from PIL import Image, ImageDraw, ImageFont
from math import *
import pandas as pd
import matplotlib.pyplot as plt



def group_read(bname,read_group):
    bamfile = pysam.AlignmentFile(bname, "rb")
    c=0
    cu=0
    for r in bamfile:
        if r.is_unmapped:
            cu+=1
            continue
        if r.is_secondary:
            continue
        if r.is_supplementary:
            continue
        if r.query_name in read_group.keys():
            read_group[r.query_name].append(r)
        else:
            read_group[r.query_name]=[r]
        c+=1
    print('mapped : {},unmapped {}'.format(str(c),str(cu)))
def compare(rg,rname,resu,FP):
    resu[rname]={'name':[],'pos':[],'ord':[]}
    max_pos=[]
    for a in rg:
        resu[rname]['name'].append(a.reference_name)
        max_pos.append(a.aend)
    result=list(set(resu[rname]['name']))
    if len(result)>2:
        print('problem')
        refpos=None
        lkey=None
    else:
        refpos,lkey=ref_loc(max_pos,result[0])
        for r in rg:
            correct=False
            resu[rname]['name']=r.reference_name
            lastpos=r.aend
            nread='peri'
            if r.tags[-1][0]=='RG': 
                resu[rname]['name']=r.reference_name
                nread=r.tags[-1][1].split('_')[1]
                correct=True
            correct_loc(r.get_reference_positions(full_length=True),nread,refpos,lastpos,FP,correct)
    return refpos,lkey


def define_loc(file,name):
    for l in file:
        tab=l.split('\t')
        if tab[0]==name:
           return int(tab[1].split(',')[1][:-1])
    return 0

def ref_loc(pos,name):
    leader=37
    lkey=[]
    file=open('data/SGtable/COV_table').readlines()[1:]
    floc=define_loc(file,name.replace('sgRNA',''))
    newpos=[]
    dic={'floc':floc,'peri':{},'COVsplice':{},'COV':{}}
    mini=0
    for j in pos:
        a=j
        if j>floc:
            a=j-(floc-37)
        if mini<a:
            mini=a
    for p in range(mini+10):
        if p<36:
            dic['peri'][str(p)]='lightgrey'
            dic['COVsplice'][str(p)]='lightgrey'
            dic['COV'][str(p)]='lightgrey'
            lkey.append(str(p))
        else:
            dic['peri'][str(p+floc-leader)]='lightgrey'
            dic['COVsplice'][str(p+floc-leader)]='lightgrey'
            dic['COV'][str(p+floc-leader)]='lightgrey'
            lkey.append(str(p+floc-leader))
    return dic,lkey


def correct_loc(pos,name,dic,lastpos,FP,correct=False):
    lastseenpos=None
    floc=0
    leader=False
    if correct:
        floc=dic['floc']-37
    for p in pos:
        if p:
            if p<36:
                dic[name][str(p)]='red'
                leader=True
            else:
                dic[name][str(p+floc)]='green'
                lastseenpos=str(p+floc)
                if leader and p<41:
                    leader=False
                    FP.append('ok')
                                
        else:
            if lastseenpos:
                dic[name][lastseenpos]='grey'
        if str(lastpos+floc)==lastseenpos:
            return
    


l=['result/minimap_COV.bam']
def main(l): 
    rg={}
    resu={}
    for f in l:
        group_read(f,rg)
    FP=[]
    for k in rg.keys():
        dic,lkey=compare(rg[k],k,resu,FP)
        # if dic:
        #     draw(dic,lkey,k)
    print(len(FP))
    


def draw(dic,lkey,name):
    width=len(dic['peri'].keys())*6
    if len(dic['peri'].keys())>6000:
        print('tropgrand')
        return
    height=600
    size = (width, height)
    image = Image.new('RGB', size,'white')
    draw = ImageDraw.Draw(image)
    font = ImageFont.truetype('Scripts/verdana.ttf', 8)
    c = 50
    px = c
    nb_rec = 1
    lh=0
    d=c+lh
    for i in ['peri','COV','COVsplice']:
        px=c
        for rec in lkey:   
            draw.rectangle([px, (height-d)-nb_rec*5, px+5, (height-d)-(nb_rec-1)*5], fill=dic[i][rec], width=1)
            px+= 5
        nb_rec += 3
    #image.show() 
    image.save('img/test_'+name+'.png')   


main(l)
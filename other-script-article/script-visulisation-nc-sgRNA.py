from BCBio import GFF
import pysam
import argparse
import pysam
from PIL import Image, ImageDraw, ImageFont
from math import *
import pandas as pd
import matplotlib.pyplot as plt
from math import ceil
pos={'ORF1ab': [265, 21555], 'S': [21562, 25384], 'ORF3a': [25392, 26220], 'E': [26244, 26472], 'M': [26522, 27191], 'ORF6': [27039, 27387], 'ORF7a': [27393, 27759], 'ORF7b': [27755, 27887], 'ORF8': [27893, 28259], 'N': [28273, 29533], 'ORF10': [29557, 29674]}
gff='covid.gff'
fm='finalresult/peri2/COV_periscope_novel_counts.csv'
fm1='all_nc_sgRNA.csv'
def main(pos,fm,fm1):
    # in_file = gff
    # in_handle = open(in_file)
    # ab= GFF.parse(in_handle)
    # end=None
    # begin=None
    # search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'
    # pos={}
    # for rec in ab:
    #     for feature in rec.features:
    #         if feature.type=='gene':
    #             begin=feature.location.nofuzzy_start
    #             end=feature.location.nofuzzy_end
    #             pos[str(feature.qualifiers['Name'][0])]=[int(begin),int(end)]
    # in_handle.close()
    # print(pos)
    # raise
    result={}
    result1={}
    for k in pos.keys():
        result[k]=[]
        result1[k]=[]
    result['other']=[]
    result1['other']=[]
    f1p=open(fm,'r').readlines()[1:]
    c=0
    cc=0
    for y in f1p:
        if y!='\n':
            ly=y.split(',')
            poso=int(ly[1].split('_')[1])
            c+=1
            for k in pos.keys():
                if pos[k][0]<poso and pos[k][1]>poso:
                    result[k].append(poso)
                    cc+=1
                    continue
            if c!=cc:
                result['other'].append(poso)
                cc+=1
    f2p=open(fm1,'r').readlines()[1:]
    print(len(f2p))
    c=0
    cc=0
    for y in f2p:
        
        if y!='\n':
            noad=False
            ly=y.split(';')
            if ly[-1]=='normal\n':
                noad=True
            if not noad:
                poso=int(ly[2])
                c+=1
                for k in pos.keys():
                    if pos[k][0]<poso and pos[k][1]>poso:
                        result1[k].append(poso)
                        cc+=1
                        break
                if c!=cc:
                    result1['other'].append(poso)
                    cc+=1
    final={}
    final1={}
    ok=[]
    for i in range(29900) :
        if i%10==0:
            ok.append(i)
            final[i]=0
            final1[i]=0

    for k in result.keys():
        for i in result[k]:
            num=int(str(i)[0:-1]+'0')
            final[num]+=1  
    for k2 in result1.keys():
        print(len(result1[k2]))
        for i in result1[k2]:
            num=int(str(i)[0:-1]+'0')
            final1[num]+=1  
    for k in ok:
        if final1[k]!=0:
            final1[k]=ceil(final1[k]/10)

        if final[k]!=0:
            final[k]=ceil(final[k]/10)

    # visualisation(final1, ok,'ok_correct.pdf',pos,height = 35000 )
    visualisation(final, ok,'ok_peri.pdf',pos,height = 35000 )
    # print(result)
    # print('n')
    # print(result1)





def visualisation(l_nucl, ok,outfile,dic,height = 55000 ):
    cov=300
    decoupe= int(sqrt(((len(l_nucl.keys())+5)*5)/((cov+5)*5)))+1
    height=int((decoupe*(cov+5)*5))-3000
    width = int((len(l_nucl.keys())*5)/decoupe)+200
    size = (width, height) 
    image = Image.new('RGB', size,'white')
    draw = ImageDraw.Draw(image)
    font = ImageFont.truetype('Scripts/verdana.ttf', 25)
    # font = ImageFont.truetype('Scripts/verdana.ttf', 8)
    c = 100
    h=1
    l=1
    px = c
    nb_rec = 1
    lh=0
    max_rec=0
    d=c+lh
    lastpx=0
    for i in ok:
        if i%500 ==0 :
            draw.text((px, height-(d-10)), f"{i}", font=font, fill ="black", align="center")
        if i%250 ==0 :
            draw.line([(px+3,height-(d-2)), (px+3,height-(d-8))], fill='black',width=5)
            
        for rec in range(l_nucl[i]+1):
            draw.rectangle([px, (height-d)-nb_rec*5, px+5, (height-d)-(nb_rec-1)*5], fill='green', width=1)
            nb_rec += 1
        if nb_rec>max_rec:
            max_rec=nb_rec
        nb_rec = 1
        px+=5
        if (i)%(int((width-50)/5)-18) == 0 and i!=0:
            ab,dic=calcline(dic,lastpx,l,c,i)
            for elem in ab:
                deb=elem[0]
                fin=elem[1]
                draw.line([(deb,height-(d-50)), (fin,height-(d-50))], fill=elem[2],width = 20)
                draw.text((deb+(fin-deb)/2, height-(d-70)), f"{elem[3]}", font=font, fill =elem[2], align="center")
                # image.show()
            lastpx=i
            px = c
            lh=5*(max_rec+5)
            d=d+c+lh
            l+=1
            
            # image.save('Comparison_Score&Coverage_Tools_to_pos{}.png'.format(i.get_pos()))
            # image.show()
            # image = Image.new('RGB', size,'white')
            # draw = ImageDraw.Draw(image)
    ab,dic=calcline(dic,lastpx,l,c,i)
    for elem in ab:
        deb=elem[0]
        fin=elem[1]
        draw.line([(deb,height-(d-50)), (fin,height-(d-50))], fill=elem[2],width = 20)
        draw.text((deb+(fin-deb)/2, height-(d-70)), f"{elem[3]}", font=font, fill =elem[2], align="center")

    image.save(outfile)
    return


def calcline(pos,wi,l,c,i):
    color={'ORF1ab': 'grey', 'S': 'blue', 'ORF3a': 'red', 'E': 'orange', 'M': 'green', 'ORF6': 'pink', 'ORF7a': 'darkblue', 'ORF7b': 'darkcyan', 'ORF8': 'brown', 'N': 'coral', 'ORF10': 'purple'}
    res=[]
    todel=[]
    ok=False
    for k in pos.keys():
        ok=False
        if pos[k][0]<i:
            if l!=1:
                if pos[k][0]==0:
                    d=c
                else:
                    d=c+int(((pos[k][0])*5)/10)-(wi/2)
            else:
                d=c+int((pos[k][0]*5)/10)
            ok=True
        if ok:
            if (pos[k][1])<=i:
                f=c+int(((pos[k][1])*5)/10)-(wi/2)
                todel.append(k)
                res.append((d,f,color[k],k))
            elif (pos[k][1])>i:
                f=c+((i-wi)/2)
                res.append((d,f,color[k],k))
                pos[k][0]=0
                break
    for j in todel:
        pos.pop(j)
    return res,pos
        
        

main(pos,fm,fm1)
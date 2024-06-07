import random

def main(linput,filename,outone,outtwo,param):
    finalfile=open(outone,'w')
    finaldata=open(outtwo,'w')
    nbread=0
    sgname=find_same_ampli(filename)
    sgpath={}
    alea=True
    for path in linput:
        f=open(path,'r').readlines()
        num=int(f[0].split('_')[0][2:])
        if num in sgname.keys():
            if sgname[num] in sgpath.keys():
                sgpath[sgname[num]].append(path)
            else:
                 sgpath[sgname[num]]=[path]
        else:
            
            if alea:
                rd=random.randint(int(param*0.8),int(param*0.9))
                nbread+=rd
            else:
                rd=param
                nbread+=rd
            for nb in range(rd):
                finalfile.write(''.join(e for e in f[nb*4:((nb*4)+4)]))
                
    lis=define_proportion(sgpath,nbread,alea)
    print(lis)
    for item in lis:
        nbitem=len(item[0])
        for path in item[0]:
            f=open(path,'r').readlines()
            finaldata.write(path+'\t'+str(int(item[1]/nbitem))+'\t')
            for nb in range(int(item[1]/nbitem)):
                finalfile.write(''.join(e for e in f[nb*4:((nb*4)+4)]))
                finaldata.write(f[nb*4].replace('\n','').replace('@','')+'\t')
            finaldata.write('\t\n')
    finalfile.close()
    finaldata.close()

def find_same_ampli(path):
    fpath=open(path).readlines()
    sgname={}
    for li in fpath:
        l=li.split('_')
        if l[0]=='>sgRNA':
            sgname[int(l[3])]=l[1]
    return sgname



def define_proportion(linput,nbread,alea):
    tpnb=list(linput.keys())
    nb=int(nbread*0.01)
    print(nb,nbread)
    tp=[]
    if alea:   
        choice=len(tpnb)-(random.randint(3,(len(tpnb)-2)))
        for _ in range(choice):
            tpnb.remove(tpnb[random.randint(0,len(tpnb)-1)])
        reste=100
        for apath in tpnb:
            choice=random.randint(int(reste*0.4),int(reste*0.75))
            reste=reste-choice
            if reste ==0:
                reste= int(nb*0.008) 
            tp.append((linput[apath],int(nb*(choice/100))))
        return tp
    else:
        choice=len(tpnb)
        for apath in tpnb:
            tp.append((linput[apath],int(100)))
        return tp



main(snakemake.input['l'],snakemake.input['p'],snakemake.output['a'],snakemake.output['b'],int(snakemake.params['nb']))
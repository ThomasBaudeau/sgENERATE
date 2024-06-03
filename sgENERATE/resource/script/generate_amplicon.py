import os
def gene_ampli_ss(ref,amplicon,out):
    ref=open(ref,'r').readlines()[1]
    fampli=open(amplicon,'r').readlines()
    f=open(out,'w')
    st_ampli=[]
    nd_ampli=[]
    ct=1
    # f.write('>nCoV-2019_0\n')
    # f.write(ref)
    for ampli in fampli:
        pos=ampli.split('\t')
        f.write('>'+pos[3])
        f.write(ref[int(pos[1]):int(pos[2])]+'\n')
        st_ampli.append(int(pos[1]))
        nd_ampli.append(int(pos[2]))
        ct+=1
    f.close()
    return st_ampli,nd_ampli,ct

def gene_ampli_sg(ref,amplicon,name,out,ct2):
    f=open(out,'a')
    for ct,ampli in enumerate(amplicon):
        f.write('>'+name.replace('sgRNA','sgRNA_')+'_'+str(ct)+'_'+str(ct2)+'\n')
        f.write(ref[ampli[0]:ampli[1]]+'\n')
        ct2+=1
    f.close()
    return ct2

def colect_pos_sg(sg,st_ampli,nd_ampli):
    fsg=open(sg,'r').readlines()
    sg=fsg[1]
    pos=fsg[0].split('_')[1].replace('\n','').strip('][').split(',')
    name=fsg[0].split('_')[0][1:]
    lst_st=[]
    new_ampli=[]
    ampsize=nd_ampli[0]-st_ampli[0]
    for ct,st in enumerate(st_ampli):
        if st<int(pos[0]):
            lst_st.append(st)
        else:
            break
    rep=binary_search(nd_ampli,int(pos[1]))
    for ct2,nd in enumerate(nd_ampli[rep:]):
        for y in lst_st:
            pos0=int(pos[0])
            pos1=int(pos[1])
            a=(pos0-y)+nd-pos1>(ampsize-int(ampsize*0.4))
            b=(int(ampsize*1.5)+int(ampsize*0.1))>(pos0-y)+(nd-pos1)
            if a and not b:
                return new_ampli,sg,name
            if a and b:
                new_ampli.append((y,nd-int(pos[1])))
    return new_ampli,sg,name



    


def binary_search(arr, x):
    low = 0
    high = len(arr) - 1
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        # If x is greater, ignore left half
        if arr[mid] < x:
            low = mid + 1
        # If x is smaller, ignore right half
        elif arr[mid] > x:
            high = mid - 1
        # means x is present at mid
        else:
            return mid
    # If we reach here, then the element was not present
    return low

def main(ss,sglist,ampli,out):
    print(ss)
    print(sglist)
    print(ampli)
    print(out)
    st_ampli,nd_ampli,ct=gene_ampli_ss(ss,ampli,out)
    for sg in sglist:
        lstNewAmpli,sg,name=colect_pos_sg(sg,st_ampli,nd_ampli)
        ct=gene_ampli_sg(sg,lstNewAmpli,name,out,ct)


main(snakemake.input['ref'],snakemake.input['c'],snakemake.input['ampli'],snakemake.output[0])
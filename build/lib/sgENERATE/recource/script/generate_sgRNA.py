import os
def make_sgRNA(ref,ltable,dos):
    out='sgRNA'
    rep=ltable.split('\t')
    os.makedirs(os.path.dirname(dos+'/'), exist_ok=True)
    with open(dos+'/'+out+rep[0]+'.fasta','w') as f:
        f.write(">"+out+rep[0]+'_'+rep[1]+'\n')
        
        stplit=int(rep[1].split(',')[0][1:])
        print(stplit)
        endplit=int(rep[1].split(',')[1][:-1])
        print(endplit)
        f.write(ref[0:stplit]+ref[endplit:])
        f.close

def main(fref,table,out):
    ref=open(fref,'r').readlines()[1]
    table=open(table,'r').readlines()[1:]
    for line in table:
        make_sgRNA(ref,line,out)

main(snakemake.input[0],snakemake.input[1],snakemake.output[0])
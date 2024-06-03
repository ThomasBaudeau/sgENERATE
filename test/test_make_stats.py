from sgENERATE.recource.script.make_stats import found_read_peri,extract_csv_info
import os

resurep={'Ground_truth': [0, 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 0], 'Periscope': [0, 98, 98, 102, 96, 100, 98, 99, 95, 0, 98, 4], 'Periscope_multi': [0, 98, 98, 93, 99, 100, 97, 99, 94, 0, 98, 1]}

path= os.path.dirname(__file__)
a=path+'/Periscope_test/'
b=path+'/Periscope_mult_test/'
c=path+'/result_test/'
file1=a+"COV_periscope_counts.csv"
file2=c+"COV_multifastq.faa"
inperi=a+"COV_periscope.bam"
filenovel=a+"COV_periscope_novel_counts.csv"
file1mult=b+"COV_periscope_counts.csv"
inperi2=b+"COV_periscope.bam"
filenovelmult=b+"COV_periscope_novel_counts.csv"
nb=c+"COV_nbread.txt"
file3=c+'final_COV_proportion.txt'
print(os.getcwd())
if file3:
    result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname =extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb,file3,LLQ=False)
else:
    result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname =extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb,LLQ=False)
found_read_peri(inperi,gpvein,list(sgcount.keys()),'peri',LLQ=False)
found_read_peri(inperi2,gpvein,list(sgcount.keys()),'peri2',LLQ=False) 

def test_nb_read():
    print(int(ttreads),36988/4)
    assert int(ttreads)==36988/4

def test_nb_sgRNA():
    TTkeys=[]
    
    for i in gpvein['GT'].keys():
        if i == 'non_canonical' or i=='ORF1a' or i=='N*':
            assert len(set(gpvein['GT'][i]))==0
        else:
            TTkeys.append(i)
            assert len(set(gpvein['GT'][i]))==100
    for i in TTkeys:
        if i == 'S':
            assert len(set(gpvein['peri'][i]))==83
            assert len(set(gpvein['peri2'][i]))==98
        if i == 'ORF3a':
            assert len(set(gpvein['peri'][i]))==95
            assert len(set(gpvein['peri2'][i]))==98
        if i == 'E':
            assert len(set(gpvein['peri'][i]))==92
            assert len(set(gpvein['peri2'][i]))==93
        if i == 'M':
            assert len(set(gpvein['peri'][i]))==87
            assert len(set(gpvein['peri2'][i]))==99
        if i == 'ORF6':
            assert len(set(gpvein['peri'][i]))==89
            assert len(set(gpvein['peri2'][i]))==100
        if i == 'ORF7a':
            assert len(set(gpvein['peri'][i]))==90
            assert len(set(gpvein['peri2'][i]))==97
        if i == 'ORF8':
            assert len(set(gpvein['peri'][i]))==93
            assert len(set(gpvein['peri2'][i]))==99
        if i == 'N':
            assert len(set(gpvein['peri'][i]))==90
            assert len(set(gpvein['peri2'][i]))==94
        if i == 'ORF10':
            assert len(set(gpvein['peri'][i]))==91
            assert len(set(gpvein['peri2'][i]))==98    
    assert len(set(gpvein['peri']['non_canonical']))==4
    assert len(set(gpvein['peri2']['non_canonical']))==1      


def test_result():
    for i in range(len(result['Periscope'])):
        assert result['Periscope'][i]==resurep['Periscope'][i]
    for i in range(len(result['Periscope_multi'])):
        assert result['Periscope_multi'][i]==resurep['Periscope_multi'][i]
    for i in range(len(result['Ground_truth'])):
        assert result['Ground_truth'][i]==resurep['Ground_truth'][i]

#test_nb_sgRNA()


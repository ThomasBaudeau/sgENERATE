from sgENERATE.recource.script.make_stats import found_read_peri,extract_csv_info
import os



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
    result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname =extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb,file3)
else:
    result,gpvein,real,ttreads,ttsgRna,ttsgRnaP,sgcount,sgRna,sgname =extract_csv_info(file1,file1mult,file2,filenovel,filenovelmult,nb)
found_read_peri(inperi,gpvein,list(sgcount.keys()),'peri')
found_read_peri(inperi2,gpvein,list(sgcount.keys()),'peri2') 



def test_nb_read():
    print(int(ttreads),36988/4)
    assert int(ttreads)==36988/4

def test_nb_sgRNA():
    for i in gpvein['GT'].keys():
        if i == 'non_canonical':
            assert len(set(gpvein['GT'][i]))==0
        else:
            assert len(set(gpvein['GT'][i]))==100
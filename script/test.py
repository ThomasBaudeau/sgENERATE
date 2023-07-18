from BCBio import GFF
a=open('data/ref/COV_ref.fna','r').readlines()[1]
b=open('COV_multifastq.faa','w')
in_file = "script/covid.gff"

in_handle = open(in_file)
ab= GFF.parse(in_handle)
end=None
begin=None
search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'
for rec in ab:
    for feature in rec.features:
        if feature.type=='gene':
            end=feature.location.nofuzzy_end
            b.write('>')
            b.write(str(feature.qualifiers['Name'][0]))
            b.write('\n')
            if begin !=None:
                if begin>feature.location.nofuzzy_start:
                    begin=feature.location.nofuzzy_start
                b.write(search+a[int(begin)-4:int(end)])
            else:
                b.write(a[0:int(end)])
            b.write('\n')
            begin=feature.location.nofuzzy_end

in_handle.close()

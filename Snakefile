NB=5000
REAL=True

def get_checkpoint_sg(wildcards):
    ck_output = checkpoints.generate_sgref.get(**wildcards).output[0]
    return expand('result/sg_ref_{species}/{sample}.fasta',species=wildcards.species, sample=glob_wildcards(os.path.join(ck_output, "{sample}.fasta")).sample)

def get_checkpoint_simu(wildcards):
    ck_output = checkpoints.simu_phase.get(**wildcards).output[0]
    return expand('result/pbsim2/{species}/simu_{nbsample}.fastq',species=wildcards.species, nbsample=glob_wildcards(os.path.join(ck_output, "simu_{nbsample}.fastq")).nbsample)


checkpoint generate_sgref:
    input:
        "data/ref/{species}_ref.fna",
        "data/SGtable/{species}_table"
    output:
        directory('result/sg_ref_{species}')  
    script:
        'script/generate_sgRNA.py'



rule generate_amplicon:
    input:
        ref="data/ref/{species}_ref.fna",
        ampli="data/amplicon/artic_amplicons_V3.bed",
        c=get_checkpoint_sg
    output:
        'result/{species}_multifastq.faa'
    script:
        'script/generate_amplicon.py'

if not REAL:
    checkpoint simu_phase:
        input:
            ref='result/{species}_multifastq.faa',
            model="data/model/R103.model"
        conda:
            "env/pbsim.yaml"
        output:
            directory("result/pbsim2/{species}/")#temp("simu_{species}")
        params:
            number=NB,
            tsize=2000,
            er='0.93',
            tmin= 2000,
            tmax= 2000
        shell:
            "mkdir {output};pbsim --depth {params.number} --prefix {output}/simu --length-mean {params.tsize} --length-min {params.tmin} --length-max {params.tmax} --accuracy-mean {params.er} --difference-ratio 23:31:46 --seed 1234 --hmm_model {input.model} {input.ref};rm {output}/*.maf;rm {output}/*.ref;"
    rule make_graph: 
        input:
            a='Periscope/{species}_periscope_counts.csv',
            peri='Periscope/{species}_periscope.bam',
            d='Periscope/{species}_periscope_novel_counts.csv',
            a2='Periscope_mult/{species}_periscope_counts.csv',
            peri2='Periscope_mult/{species}_periscope.bam',
            d2='Periscope_mult/{species}_periscope_novel_counts.csv',
            b='result/{species}_multifastq.faa',
            c='result/final_{species}_proportion.txt',
            nbread='result/{species}_nbread.txt'
        output:
            'result/final_{species}_SG_ER.pdf'
        conda:
            'env/pyt.yaml'
        shell:
            'touch {output}'
        # script:
        #     'script/make_stats.py'

    rule agregate_read:
        input:
            l=get_checkpoint_simu,
            p='result/{species}_multifastq.faa',
        output:
            a='result/final_{species}_agregate.fastq',
            b='result/final_{species}_proportion.txt'
        params:
            nb=NB
        script:
            'script/agregate_amplicon.py'
else:
    rule make_graph: 
        input:
            a='Periscope/{species}_periscope_counts.csv',
            peri='Periscope/{species}_periscope.bam',
            d='Periscope/{species}_periscope_novel_counts.csv',
            a2='Periscope_mult/{species}_periscope_counts.csv',
            peri2='Periscope_mult/{species}_periscope.bam',
            d2='Periscope_mult/{species}_periscope_novel_counts.csv',
            b='result/{species}_multifastq.faa',
            nbread='result/{species}_nbread.txt'
        output:
            'result/final_{species}_SG_ER.pdf'
        conda:
            'env/pyt.yaml'
        shell:
            'touch {output}'
        # script:
        #     'script/make_stats.py'


rule periscope:
    input:
        'result/final_{species}_agregate.fastq'
    output:
        'Periscope/{species}_periscope_counts.csv',
        'Periscope/{species}_periscope_novel_counts.csv',
        'Periscope/{species}_periscope.bam'
    conda:
        'env/peri.yaml'
    benchmark:
        "benchmarks/Periscope_{species}.txt"
    shell:
        """
        periscope --fastq {input} --technology ont --artic-primers V3 --output-prefix {wildcards.species} --threads 8;
        mv {wildcards.species}_periscope_counts.csv Periscope/{wildcards.species}_periscope_counts.csv;
        mv {wildcards.species}_periscope_novel_counts.csv Periscope/{wildcards.species}_periscope_novel_counts.csv
        mv {wildcards.species}_periscope.bam Periscope/{wildcards.species}_periscope.bam
        rm COV*
        """ 

rule all:
    input:
        'result/final_COV_SG_ER.pdf'


rule run_minimap2:
    input:
        ref='result/{species}_multifastq.faa',
        fastq='result/final_{species}_agregate.fastq'
    output:
        'result/minimap_{species}.bam'
    conda:
        'env/minimap2.yaml'
    shell:
        "minimap2 -ax map-ont --secondary=no -k 15  -w 1 -t 8 {input.ref} {input.fastq} | samtools sort - | samtools view -bh  > {output};samtools index {output}"
        #minimap2 -ax splice -k 8 -w 1 -g 20000 -G 20000 -CO -uf --no-end-flt --splice-flank=no {input.ref} {input.fastq} | samtools sort - | samtools view -bh  > {output};samtools index {output}"

rule run_periscope_multi:
    input:
        'result/final_{species}_agregate.fastq'
    output:
        'Periscope_mult/{species}_periscope_counts.csv',
        'Periscope_mult/{species}_periscope_novel_counts.csv',
        'Periscope_mult/{species}_periscope.bam'
    conda:
        'env/perimult.yaml'
    benchmark:
        "benchmarks/PeriscopeMult_{species}.txt
    shell:
        "periscope --fastq {input} --gff script/covid.gff --technology ont --artic-primers V3 --output-prefix Periscope_mult/{wildcards.species} --threads 8;rm Periscope_mult/COV.*"

rule getnbread:
    input:
        'result/final_{species}_agregate.fastq'
    output:
        'result/{species}_nbread.txt'
    shell:
        'cat {input} | wc -l > {output}'
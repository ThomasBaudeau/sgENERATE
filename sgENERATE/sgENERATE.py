import argparse
import sys
import os
import snakemake
import glob
import logging

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='sgENERATE: a tool to create ARTIC like nanopore data and compare sgRNA finding tools',usage='''sgENERATE [options]''')
    parser.add_argument('--coverage',dest='cov', help='Coverage of the generated samples', default=5000,required=False)
    parser.add_argument('--error',dest='er' ,help='Mean error rate of the generated samples (default 0.95, i.e., 5 percent error)',default=0.95,required=False)
    parser.add_argument('--real',dest='real',help='Location of the fastq file to be used for tools benchmarking',required=False,default=False)
    parser.add_argument('--no-compare',dest='comp', help='Default False, use it if only dataset generation is desired',action='store_false',required=False)
    parser.add_argument('--tool',dest='tool', help='Advanced option: see online documentation for more information',default='minimap',required=False)
    #parser.add_argument('--artic-primers', dest='artic_primers', help='artic network primer version used:\n* V1 (default), V2, V3, V4\n* 2kb (for the UCL longer amplicons)\n* midnight (1.2kb midnight amplicons)\n* for custom primers provide path to amplicons file first and primers file second', nargs='*', default="V3")
    parser.add_argument('--technology', help='the sequencing technology used, either:\n*ont\n*illumina', default="ont")
    parser.add_argument('--fastq',dest='fastq',help='if you already have a single fastq then you can use this flag instead, if illumina paired end separate fastq by space', nargs='+',required=False,default=['result/final_{species}_agregate.fastq'])
    parser.add_argument('--mode',dest='mode',help='use --mode extraction for extract all the non-canonical sgRNA',required=False,default=False)
    parser.set_defaults(comp=True)
    args = parser.parse_args()


    # run snakemake pipeline 1st
    dir = os.path.join(os.path.dirname(__file__))
    scripts_dir= os.path.join(dir, 'resource')

    config = dict(
        NB=args.cov,
        ER=args.er,
        REAL=args.real,
        COMP=args.comp,
        tool=args.tool,
        path=scripts_dir,
        tech=args.technology,
        fastq=args.fastq,
        mode=args.mode
    )

    print(config)

    snakefile = os.path.join(scripts_dir, 'Snakefile')
    print(snakefile)
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)
    else:
        print("Found the snakefile")

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=False, forceall=True, force_incomplete=True,
                                 config=config, cores=16, lock=False, use_conda=True,
                                 )
    if status:  # translate "success" into shell exit code of 0
        exit(0)
    else: 
        'error'
        exit(1)

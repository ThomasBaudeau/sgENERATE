import argparse
import sys
import os
import snakemake
import glob
import logging

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='sgENERATE: a tool to create ARTIC like nanopore data and compare sgRNA finding tools',usage='''periscope [options]''')
    parser.add_argument('--coverage',dest='cov', help='coverage of the generated samples', default=5000,required=False)
    parser.add_argument('--error',dest='er' ,help='error rate of the samples',default=0.95,required=False)
    parser.add_argument('--real',dest='real',help='provide a real dataset for benching the tools',required=False,default=False)
    parser.add_argument('--compare',dest='comp', help='make also the comparaison',default=True,required=False)
    #parser.add_argument('--artic-primers', dest='artic_primers', help='artic network primer version used:\n* V1 (default), V2, V3, V4\n* 2kb (for the UCL longer amplicons)\n* midnight (1.2kb midnight amplicons)\n* for custom primers provide path to amplicons file first and primers file second', nargs='*', default="V3")



    args = parser.parse_args()


    # run snakemake pipeline 1st
    dir = os.path.join(os.path.dirname(__file__))
    scripts_dir= dir

    config = dict(
        NB=args.cov,
        ER=args.er,
        REAL=args.real,
        COMP=args.comp
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

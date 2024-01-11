#!/usr/bin/env python3
from Bio import pairwise2
import pysam
import argparse
from pybedtools import *
from artic.vcftagprimersites import read_bed_file
import sys
from numpy import median
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor as ProcessPool, process
import time

class ClassifiedRead():
    def __init__(self,sgRNA: bool,orf: str,read: pysam.AlignedRead):
        self.sgRNA = sgRNA
        self.orf = orf
        self.pos = read.pos
        self.read = read.to_string()

def get_mapped_reads(bam):
    mapped_reads = int(pysam.flagstat(bam).split("\n")[4].split(" ")[0])-int(pysam.flagstat(bam).split("\n")[2].split(" ")[0])-int(pysam.flagstat(bam).split("\n")[1].split(" ")[0])
    return mapped_reads

# def check_start(bed_object,read):
#     """
#     find out where the read is in a bed file, in this case the ORF starts
#     :param bed_object: bedtools object
#     :param read: pysam read object
#     :return: the orf
#     """

#     # reads with a pos of 0 make this fail so puting in try except works
#     try:
#         read_feature = BedTool(read.reference_name + "\t" + str(read.pos) + "\t" + str(read.pos), from_string=True)
#         intersect = bed_object.intersect(read_feature)
#         orf=intersect[0].name
#         # if len(intersect) > 1:
#         #     print("odd")
#     except:
#         orf=None
#     # remove bedtools objects from temp
#     cleanup()
#     return orf

def check_start(read, leader_search_result, orfBed):
    orf=None
    for row in orfBed:
        # see if read falls within ORF start location
        if row.end >= read.reference_start >= row.start:
            orf = row.name
    if orf == None:
        if leader_search_result == True:
            orf = "novel_" + str(read.reference_start)
    return orf




def extact_soft_clipped_bases(read):
    """
            +-----+--------------+-----+
            |M    |BAM_CMATCH    |0    |
            +-----+--------------+-----+
            |I    |BAM_CINS      |1    |
            +-----+--------------+-----+
            |D    |BAM_CDEL      |2    |
            +-----+--------------+-----+
            |N    |BAM_CREF_SKIP |3    |
            +-----+--------------+-----+
            |S    |BAM_CSOFT_CLIP|4    |
            +-----+--------------+-----+
            |H    |BAM_CHARD_CLIP|5    |
            +-----+--------------+-----+
            |P    |BAM_CPAD      |6    |
            +-----+--------------+-----+
            |=    |BAM_CEQUAL    |7    |
            +-----+--------------+-----+
            |X    |BAM_CDIFF     |8    |
            +-----+--------------+-----+
            |B    |BAM_CBACK     |9    |
            +-----+--------------+-----+
    :param cigar:
    :return:
    """


    cigar = read.cigartuples
    if cigar[0][0] == 4:
        # there is softclipping on 5' end
        number_of_bases_sclipped = cigar[0][1]
        # get 3 extra bases incase of homology
        bases_sclipped = read.seq[0:number_of_bases_sclipped+3]



        search='AACCAACTTTCGATCTCTTGTAGATCTGTTCTC'


        # not enough bases to determine leader
        if number_of_bases_sclipped < 6:
            return False

        # determine perfect score for a match
        if number_of_bases_sclipped >= 33:
            # allow 3 mistamches
            perfect = 60.0
        else:
            # allow 1 mistamches
            perfect = (number_of_bases_sclipped * 2)-2

        align = pairwise2.align.localms(bases_sclipped, search, 2, -2, -20, -.1,  one_alignment_only=True  )
        align_score=align[0][2]
        align_right_position = align[0][4]

        logger.debug("Processing read: %s", read.query_name)
        logger.debug("%s perfect score: %s", read.query_name,str(perfect))
        logger.debug("%s actual score: %s", read.query_name, str(align_score))
        logger.debug("%s alignment: %s", read.query_name, align)

        # position of alignment must be all the way to the right
        if align_right_position >= len(search):
            # allow only one mismtach - Mismatches already taken care of above?
            # if perfect-align_score <= 2:
            if perfect-align_score <= 0:
                logger.debug("%s sgRNA: %s", read.query_name, True)
                return True
            else:
                # more than allowed number of mismatches
                logger.debug("%s sgRNA: %s,%s", read.query_name, False,'too many mismatches')
                return False
        else:
            # match is not at the end of the leader
            logger.debug("%s sgRNA: %s,%s", read.query_name, False,'match not at end of leader')
            return False
    else:
        # No soft clipping at 5' end of read
        logger.debug("%s sgRNA: %s,%s", read.query_name, False,'no soft-clipping')
        return False






def open_bed(bed):
    """
    open bed file and return a bedtools object
    :param bed:
    :return:
    """
    bed_object = BedTool(bed)
    return bed_object



def process_reads(data):
    bam = data[0]
    args = data[1]
    inbamfile = pysam.AlignmentFile(bam, "rb")
    bam_header = inbamfile.header.copy().to_dict()
    outbamfile = pysam.AlignmentFile(bam + "_periscope_temp.bam", "wb", header=bam_header)
    read_dic={}
    mapped_reads = get_mapped_reads(bam)
    logger.warning("Processing " + str(mapped_reads) + " reads")

    orf_bed_object = open_bed(args.orf_bed)

    reads={}
    for read in inbamfile:

        if read.seq == None:
            # print("%s read has no sequence" %
            #       (read.query_name), file=sys.stderr)
            continue
        if read.is_unmapped:
            # print("%s skipped as unmapped" %
            #       (read.query_name), file=sys.stderr)
            continue
        if read.is_supplementary:
            # print("%s skipped as supplementary" %
            #       (read.query_name), file=sys.stderr)
            continue  
        if read.is_secondary:
            # print("%s skipped as secondary" %
            #       (read.query_name), file=sys.stderr)
            continue
        if read.query_name in read_dic:
            continue
        # print("------")
        # print(read.query_name)
        # # print(read.is_read1)
        # # print(read.get_tags())
        # print(read.cigar)
        leader_search_result = extact_soft_clipped_bases(read)

        if read.query_name not in reads:
            reads[read.query_name] = []

        orfRead = check_start(read, leader_search_result, orf_bed_object)
        reads[read.query_name].append(

            ClassifiedRead(sgRNA=leader_search_result,orf=orfRead,read=read)


        )
        read.set_tag('XS', '/')
        read.set_tag('XA', '/')
        if leader_search_result:
            read_dic[read.query_name]=0
            if 'novel_' in orfRead:
                read.set_tag('XC', orfRead)
            else:
                read.set_tag('XC', 'sgRNA')
                read.set_tag('XO', orfRead)
        else:
            read.set_tag('XC', 'gRNA')



        # ok now add this info to a dictionary for later processing


        # write the annotated read to a bam file
        outbamfile.write(read)
    outbamfile.close()    
    return(reads)

def process_pairs(reads_dict):
    # now we have all the reads classified, deal with pairs
    logger.info("dealing with read pairs")

    orfs={}
    orfs_gRNA={}
    for id,pair in reads_dict.items():
        # get the class and orf of the left hand read, this will be the classification and ORF for the pair - sometimes right read looks like it has subgenomic evidence - there are likely false positives

        left_read = min(pair, key=lambda x: x.pos)

        read_class = left_read.sgRNA
        orf = left_read.orf

        left_read_object = left_read.read

        if orf == None:
            continue

        #build orfs_gRNA for every non-novel sgRNA
        if orf not in orfs_gRNA and "novel" not in orf:
            orfs_gRNA[orf] = 0

        #assign read to gRNA
        if read_class == False:
            orfs_gRNA[orf] += 1
        else:
            #assign read to sgRNA
            if orf not in orfs:
                orfs[orf] = [left_read_object]
            else:
                orfs[orf].append(left_read_object)

    logger.info("dealing with read pairs....DONE")

    return orfs, orfs_gRNA

def multiprocessing(func, args, workers):
    with ProcessPool(workers) as ex:
        res = list(tqdm(ex.map(func, args), total=len(args)))
    return res

def main(args):

    # t1=time.time()
    inbamfile = pysam.AlignmentFile(args.bam, "rb")
    # bam_header = inbamfile.header.copy().to_dict()

    #get a list of bams:
    import glob
    files = glob.glob(args.output_prefix+".split.*.sam")
    threads=args.threads
    result=[]
    for file in files:
        result.append([file,args])

    processed = multiprocessing(
        process_reads,
        args=result,
        workers=int(threads)
    )
    # outbamfile.close()
    output_bams = [file+"_periscope_temp.bam" for file in files]
    pysam.merge(*["-f",args.output_prefix + "_periscope.bam"]+output_bams)
    # pysam.sort("-o", args.output_prefix + "_periscope_sorted.bam",  args.output_prefix + "_periscope.bam")
    # pysam.index(args.output_prefix + "_periscope_sorted.bam")

    # t2=time.time()
    # print("periscope.py time:", t2-t1)
# print(os.getcwd())
# process_reads(['output_periscope.bam','sgENERATE/recource/data/SGtable/orf_start.bed'])
if __name__ == '__main__':



    parser = argparse.ArgumentParser(description='periscopre: Search for sgRNA reads in artic network SARS-CoV-2 sequencing data')
    parser.add_argument('--bam', help='bam file',default="The bam file of full artic reads")
    parser.add_argument('--output-prefix',dest='output_prefix',help="Path to the output, e.g. <DIR>/<SAMPLE_NAME>")
    parser.add_argument('--threads',dest='threads',help="Path to the output, e.g. <DIR>/<SAMPLE_NAME>")
    parser.add_argument('--orf-bed', dest='orf_bed', help='The bed file with ORF start positions')

    logger = logging
    logger.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    logger.info("staring periscope")

    args = parser.parse_args()

    periscope = main(args)

    if periscope:
        print("all done", file=sys.stderr)





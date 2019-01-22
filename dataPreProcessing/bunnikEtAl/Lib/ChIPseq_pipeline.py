import sys
import tempfile
import os
import random

#The first 5 bases and the last base were systematically removed from the
#sequence reads using FastQ Trimmer, part of the FASTX-Toolkit
#(http://hannonlab.cshl.edu/fastx_toolkit/index.html). Contaminating adaptor
#reads were removed using Scythe
#(https://github.com/ucdavis-bioinformatics/scythe). Reads were then trimmed for
#bases with a quality score below 30, and reads containing any Ns as well as
#reads shorter than 18 bases were discarded using Sickle
#(https://github.com/ucdavis-bioinformatics/sickle). The trimmed sequence reads
#were first mapped to the human genome (HG19, downloaded from
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/), and all non-mapped reads were
#subsequently mapped to P. falciparum 3D7 genome v9.0 (downloaded from
#http://www.plasmoDB.org) using BWA [52] with default error rates, allowing a
#maximum of 1500 bp distance between read pairs. Any read that was either
#non-uniquely mapped (Samtools v0.1.18 [53]), not properly paired (Samtools) or
#a PCR duplicate (Picard Tools v1.78 [http://picard.sourceforge.net/]) was
#discarded. The final number of mapped reads for each library is listed in
#Additional file 1: Table S1.

adapter_file = "/net/noble/vol2/home/katecook/software/scythe-master"


def main():
    #random.seed(1)
    
    fastq_file = sys.argv[1]
    #trim_length = int(sys.argv[2])
    bwa_index = sys.argv[2]
    
    basename = os.path.basename(fastq_file)
    prefix = (basename.split("."))[0]
    
    # trim first 5 bases 
    #tempdir = tempfile.gettempdir()
    tempdir = "/net/noble/vol2/home/katecook/proj/2016predictExpression/results/katecook/20161119_bunnik2014_MNase_data"
    trimmed_file =  "%s/%s.%d" % (tempdir, basename, random.randint(1,1000))
    trim_cmd = "zcat %s | fastx_trimmer -f 6 -z -o %s" % (fastq_file,trimmed_file)
    print trim_cmd
    os.system(trim_cmd)
    
    # remove adapters
    scythed_file = "%s_scythed" % trimmed_file
    scythe_cmd = "scythe -a %s %s > %s" % (adapter_file, trimmed_file, scythed_file)
    os.system(trim_cmd)
    
    sys.exit()
    
    # filter reads
    sickle_cmd = "sickle pe -f" 
    
    # map sequences to supplied genome
    sam_file = "%s/%s.sam" % (tempdir, prefix)
    map_cmd = "bwa mem %s %s > %s" % (bwa_index,trimmed_file,sam_file)
    print map_cmd
    os.system(map_cmd)

    # convert to bam & sort
    bam_file = "%s/%s.bam" % (tempdir, prefix)
    sort_tempfile_prefix = "%s/%s.%d" % (tempdir, prefix, random.randint(1,1000))
    bam_cmd = "samtools sort -o %s -O bam -T %s %s" % (bam_file, sort_tempfile_prefix, sam_file)
    print bam_cmd
    os.system(bam_cmd)
    
    # index bamfile
    index_cmd = "samtools index %s" % bam_file
    index_file = bam_file + ".bai"
    print index_cmd
    os.system(index_cmd)
    
    # filter mapped reads
    filtered_bam_file = "%s_filtered.bam" % prefix
    filter_cmd =  "samtools view -F 1804 -b -o %s %s" % (filtered_bam_file, bam_file)
    print filter_cmd
    os.system(filter_cmd)
    
    # delete temporary files
    os.remove(trimmed_file)
    os.remove(sam_file)
    os.remove(bam_file)
    os.remove(index_file)

if __name__ == "__main__":
    main()

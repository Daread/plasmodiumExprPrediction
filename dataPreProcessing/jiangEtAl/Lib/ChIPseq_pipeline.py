import sys
import tempfile
import os
import random

def main():
    #random.seed(1)
    
    fastq_file = sys.argv[1]
    trim_length = int(sys.argv[2])
    bwa_index = sys.argv[3]
    
    basename = os.path.basename(fastq_file)
    prefix = (basename.split("."))[0]
    
    # trim sequences to supplied length
    tempdir = tempfile.gettempdir()
    #tempdir = "/net/noble/vol2/home/katecook/proj/2016predictExpression/results/katecook/20161115_bartfai2010_histone_data"
    trimmed_file =  "%s/%s.%d" % (tempdir, basename, random.randint(1,1000))
    trim_cmd = "zcat %s | fastx_trimmer -l %d -z -o %s" % (fastq_file,trim_length,trimmed_file)
    print trim_cmd
    os.system(trim_cmd)
    
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

FILES = Nucleosome_0h.bam  Nucleosome_18h.bam  Nucleosome_36h.bam

BWA_INDEX = /net/noble/vol1/data/bwa-indices/pfal3D7-PlasmoDBv29/pfal3D7-PlasmoDBv29

echo:
	echo $(ACCESSIONS)
	echo $(addsuffix .fastq.gz, $(ACCESSIONS))

download_fastq:
	$(foreach a, $(ACCESSIONS), \
	   fastq-dump --gzip $(a) $(a).fastq.gz; \
	)

fastqc:
	fastqc $(addsuffix .fastq.gz, $(ACCESSIONS))
	mkdir -p fastqc
	mv *_fastqc.html *_fastqc.zip fastqc

test_pipeline:
	python Lib/ChIPseq_pipeline.py SRR923450.fastq.gz $(BWA_INDEX)

run_pipeline:
	$(foreach a, $(ACCESSIONS), \
	   echo "#!/bin/bash" > $(a).job; \
	   echo "#$$ -cwd" >> $(a).job; \
	   echo "#$$ -V" >> $(a).job; \
	   echo "#$$ -l h_vmem=8G" >> $(a).job; \
	   echo "#$$ -l h_rt=36000" >> $(a).job; \
	   echo "source ~/myenv.sh" >> $(a).job; \
	   echo hostname >> $(a).job; \
	   echo "python Lib/ChIPseq_pipeline.py $(a).fastq.gz $(SEQ_TRIM_LENGTH) $(BWA_INDEX)" >> $(a).job; \
	   qsub $(a).job; \
	)

delete_fastq:
	$(foreach a, $(ACCESSIONS), \
	   rm -f $(a).fastq.gz; \
	)

create_bedfiles:
	$(foreach a, $(FILES), \
	   echo "#!/bin/bash" > $(a).job; \
	   echo "#$$ -cwd" >> $(a).job; \
	   echo "#$$ -V" >> $(a).job; \
	   echo "#$$ -l h_vmem=8G" >> $(a).job; \
	   echo "#$$ -l h_rt=36000" >> $(a).job; \
	   echo "source ~/myenv.sh" >> $(a).job; \
	   echo hostname >> $(a).job; \
	   echo "genomeCoverageBed -bg -ibam $(a) > $(a).bedgraph.tmp" >> $(a).job; \
	   echo "cat $(a).bedgraph.tmp | sed 's/GeneDB|//g' > $(a).bedgraph" >> $(a).job; \
	   echo "rm -f  $(a).bedgraph.tmp" >> $(a).job; \
	   qsub $(a).job; \
	)

IDS = Nucleosome_0h Nucleosome_18h  Nucleosome_36h

split_by_insert_length:
	$(foreach a, $(IDS), \
	   echo "#!/bin/bash" > $(a).job; \
	   echo "#$$ -cwd" >> $(a).job; \
	   echo "#$$ -V" >> $(a).job; \
	   echo "#$$ -l h_vmem=8G" >> $(a).job; \
	   echo "#$$ -l h_rt=36000" >> $(a).job; \
	   echo "source ~/myenv.sh" >> $(a).job; \
	   echo hostname >> $(a).job; \
	   echo "samtools view $(a).bam > $(a).sam;" >> $(a).job; \
	   echo "samtools view -H $(a).bam > $(a)_header.sam;" >> $(a).job; \
	   echo "python Lib/split_sam_by_insert_length.py $(a).sam $(a);" >> $(a).job; \
	   echo "rm -f $(a).sam;" >> $(a).job; \
	   echo "cat $(a)_header.sam $(a)_lt_100.sam | samtools view -Sb - > $(a)_lt_100.bam;" >> $(a).job; \
	   echo "cat $(a)_header.sam $(a)_100_to_200.sam | samtools view -Sb - > $(a)_100_to_200.bam;" >> $(a).job; \
	   echo "cat $(a)_header.sam $(a)_gt_200.sam | samtools view -Sb - > $(a)_gt_200.bam;" >> $(a).job; \
	   echo "rm -f $(a)_header.sam $(a)_lt_100.sam $(a)_100_to_200.sam $(a)_gt_200.sam;" >> $(a).job; \
	   echo "genomeCoverageBed -bg -ibam $(a)_lt_100.bam > $(a)_lt_100.bedgraph" >> $(a).job; \
	   echo "genomeCoverageBed -bg -ibam $(a)_100_to_200.bam > $(a)_100_to_200.bedgraph" >> $(a).job; \
	   echo "genomeCoverageBed -bg -ibam $(a)_gt_200.bam > $(a)_gt_200.bedgraph" >> $(a).job; \
	   qsub $(a).job; \
	)

macs:
	echo "#!/bin/bash" > header.job; \
	echo "#$$ -cwd" >> header.job; \
	echo "#$$ -V" >> header.job; \
	echo "#$$ -l h_vmem=8G" >> header.job; \
	echo "#$$ -l h_rt=36000" >> header.job; \
	echo "source ~/myenv.sh" >> header.job; \
	echo hostname >> header.job;
	cp header.job .job;



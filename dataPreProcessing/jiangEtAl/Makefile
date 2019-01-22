ACCESSIONS = $(shell cat accession_list.txt)

SEQ_TRIM_LENGTH = 50 

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
	python Lib/ChIPseq_pipeline.py SRR065648.fastq.gz $(SEQ_TRIM_LENGTH) $(BWA_INDEX)

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
	$(foreach a, $(ACCESSIONS), \
	   echo "#!/bin/bash" > $(a).job; \
	   echo "#$$ -cwd" >> $(a).job; \
	   echo "#$$ -V" >> $(a).job; \
	   echo "#$$ -l h_vmem=8G" >> $(a).job; \
	   echo "#$$ -l h_rt=36000" >> $(a).job; \
	   echo "source ~/myenv.sh" >> $(a).job; \
	   echo hostname >> $(a).job; \
	   echo "genomeCoverageBed -bg -ibam $(a)_filtered.bam > $(a).bedgraph" >> $(a).job; \
	   qsub $(a).job; \
	)




ACCESSIONS = $(shell cat accession_list.txt)
ACCESSIONS_MERGED = $(shell cat accession_list.txt | grep -v SRR065662 | grep -v SRR065663 | cap.pl SRR065662_3 )

SEQ_TRIM_LENGTH = 76

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

test_wiggler:
	align2rawsignal -i=SRR065647_filtered.bam -s=/net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/chromosomes_Pf

merge:
	samtools merge SRR065662_3_filtered.bam SRR065662_filtered.bam SRR065663_filtered.bam

create_bedfiles:
	$(foreach a, $(ACCESSIONS_MERGED), \
	   echo "#!/bin/bash" > $(a).job; \
	   echo "#$$ -cwd" >> $(a).job; \
	   echo "#$$ -V" >> $(a).job; \
	   echo "#$$ -l h_vmem=8G" >> $(a).job; \
	   echo "#$$ -l h_rt=36000" >> $(a).job; \
	   echo "source ~/myenv.sh" >> $(a).job; \
	   echo hostname >> $(a).job; \
	   echo "samtools sort $(a)_filtered.bam -o $(a)_sorted.bam -T $(a)" >> $(a).job; \
	   echo "genomeCoverageBed -ibam $(a)_sorted.bam -g > $(a).bed" >> $(a).job; \
	   echo "genomeCoverageBed -bg -ibam $(a)_sorted.bam > $(a).bedgraph" >> $(a).job; \
	   qsub $(a).job; \
	)

index_bams:
	$(foreach a, $(ACCESSIONS_MERGED), \
	   echo "#!/bin/bash" > $(a).job; \
	   echo "#$$ -cwd" >> $(a).job; \
	   echo "#$$ -V" >> $(a).job; \
	   echo "#$$ -l h_vmem=8G" >> $(a).job; \
	   echo "#$$ -l h_rt=36000" >> $(a).job; \
	   echo "source ~/myenv.sh" >> $(a).job; \
	   echo hostname >> $(a).job; \
	   echo "samtools index $(a)_sorted.bam " >> $(a).job; \
	   qsub $(a).job; \
	)

foldchange_over_input_H2A.Z:
	echo "#!/bin/bash" > foldchange_H2A.Z.job;
	echo "#$$ -cwd" >> foldchange_H2A.Z.job;
	echo "#$$ -V" >> foldchange_H2A.Z.job;
	echo "#$$ -l h_vmem=8G" >> foldchange_H2A.Z.job;
	echo "#$$ -l h_rt=36000" >> foldchange_H2A.Z.job;
	echo "source ~/myenv.sh" >> foldchange_H2A.Z.job;
	echo hostname >> foldchange_H2A.Z.job;
	echo "bamCompare -b1 SRR065647_sorted.bam -b2 SRR065660_sorted.bam -o ChIPseq_H2A.Z_10hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H2A.Z.job;
	echo "bamCompare -b1 SRR065648_sorted.bam -b2 SRR065661_sorted.bam -o ChIPseq_H2A.Z_20hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H2A.Z.job;
	echo "bamCompare -b1 SRR065649_sorted.bam -b2 SRR065662_3_sorted.bam -o ChIPseq_H2A.Z_30hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H2A.Z.job;
	echo "bamCompare -b1 SRR065650_sorted.bam -b2 SRR065664_sorted.bam -o ChIPseq_H2A.Z_40hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H2A.Z.job;
	qsub foldchange_H2A.Z.job;

foldchange_over_input_H3K9ac:
	echo "#!/bin/bash" > foldchange_H3K9ac.job;
	echo "#$$ -cwd" >> foldchange_H3K9ac.job;
	echo "#$$ -V" >> foldchange_H3K9ac.job;
	echo "#$$ -l h_vmem=8G" >> foldchange_H3K9ac.job;
	echo "#$$ -l h_rt=36000" >> foldchange_H3K9ac.job;
	echo "source ~/myenv.sh" >> foldchange_H3K9ac.job;
	echo hostname >> foldchange_H3K9ac.job;
	echo "bamCompare -b1 SRR065652_sorted.bam -b2 SRR065660_sorted.bam -o ChIPseq_H3K9ac_10hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K9ac.job;
	echo "bamCompare -b1 SRR065653_sorted.bam -b2 SRR065661_sorted.bam -o ChIPseq_H3K9ac_20hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K9ac.job;
	echo "bamCompare -b1 SRR065654_sorted.bam -b2 SRR065662_3_sorted.bam -o ChIPseq_H3K9ac_30hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K9ac.job;
	echo "bamCompare -b1 SRR065655_sorted.bam -b2 SRR065664_sorted.bam -o ChIPseq_H3K9ac_40hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K9ac.job;
	qsub foldchange_H3K9ac.job;

foldchange_over_input_H3K4me3:
	echo "#!/bin/bash" > foldchange_H3K4me3.job;
	echo "#$$ -cwd" >> foldchange_H3K4me3.job;
	echo "#$$ -V" >> foldchange_H3K4me3.job;
	echo "#$$ -l h_vmem=8G" >> foldchange_H3K4me3.job;
	echo "#$$ -l h_rt=36000" >> foldchange_H3K4me3.job;
	echo "source ~/myenv.sh" >> foldchange_H3K4me3.job;
	echo hostname >> foldchange_H3K4me3.job;
	echo "bamCompare -b1 SRR065656_sorted.bam -b2 SRR065660_sorted.bam -o ChIPseq_H3K4me3_10hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K4me3.job;
	echo "bamCompare -b1 SRR065657_sorted.bam -b2 SRR065661_sorted.bam -o ChIPseq_H3K4me3_20hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K4me3.job;
	echo "bamCompare -b1 SRR065658_sorted.bam -b2 SRR065662_3_sorted.bam -o ChIPseq_H3K4me3_30hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K4me3.job;
	echo "bamCompare -b1 SRR065659_sorted.bam -b2 SRR065664_sorted.bam -o ChIPseq_H3K4me3_40hpi_log2ratio.bedgraph -of bedgraph --scaleFactorsMethod readCount --ratio log2 --binSize 100" >> foldchange_H3K4me3.job;
	qsub foldchange_H3K4me3.job;



#macs:
#	echo "#!/bin/bash" > header.job; \
#    echo "#$$ -cwd" >> header.job; \
#    echo "#$$ -V" >> header.job; \
#    echo "#$$ -l h_vmem=8G" >> header.job; \
#    echo "#$$ -l h_rt=36000" >> header.job; \
#    echo "source ~/myenv.sh" >> header.job; \
#    echo hostname >> header.job;
#	cp header.job SRR065647.job;
#	echo "macs2 bdgcmp -t SRR065647.bedgraph -c SRR065660.bedgraph -o SRR065647_fc.bedgraph.tmp -m subtract" >> SRR065647.job;
#	echo "cat SRR065647_fc.bedgraph.tmp | grep -v nan > SRR065647_fc.bedgraph" >> SRR065647.job;
#	qsub SRR065647.job;
#	cp header.job SRR065648.job;
#	echo "macs2 bdgcmp -t SRR065648.bedgraph -c SRR065661.bedgraph -o SRR065648_fc.bedgraph.tmp -m subtract" >> SRR065648.job;
#	echo "cat SRR065648_fc.bedgraph.tmp | grep -v nan > SRR065648_fc.bedgraph" >> SRR065648.job;
#	qsub SRR065648.job;
#	cp header.job SRR065649.job;
#	echo "macs2 bdgcmp -t SRR065649.bedgraph -c SRR065662_3.bedgraph -o SRR065649_fc.bedgraph.tmp -m subtract" >> SRR065649.job;
#	echo "cat SRR065649_fc.bedgraph.tmp | grep -v nan > SRR065649_fc.bedgraph" >> SRR065649.job;
#	qsub SRR065649.job;
#	cp header.job SRR065650.job;
#	echo "macs2 bdgcmp -t SRR065650.bedgraph -c SRR065664.bedgraph -o SRR065650_fc.bedgraph.tmp -m subtract" >> SRR065650.job;
#	echo "cat SRR065650_fc.bedgraph.tmp | grep -v nan > SRR065650_fc.bedgraph" >> SRR065650.job;
#	qsub SRR065650.job;
#	cp header.job SRR065652.job;
#	echo "macs2 bdgcmp -t SRR065652.bedgraph -c SRR065660.bedgraph -o SRR065652_fc.bedgraph.tmp -m subtract" >> SRR065652.job;
#	echo "cat SRR065652_fc.bedgraph.tmp | grep -v nan > SRR065652_fc.bedgraph" >> SRR065652.job;
#	qsub SRR065652.job;
#	cp header.job SRR065653.job;
#	echo "macs2 bdgcmp -t SRR065653.bedgraph -c SRR065661.bedgraph -o SRR065653_fc.bedgraph.tmp -m subtract" >> SRR065653.job;
#	echo "cat SRR065653_fc.bedgraph.tmp | grep -v nan > SRR065653_fc.bedgraph" >> SRR065653.job;
#	qsub SRR065653.job;
#	cp header.job SRR065654.job;
#	echo "macs2 bdgcmp -t SRR065654.bedgraph -c SRR065662_3.bedgraph -o SRR065654_fc.bedgraph.tmp -m subtract" >> SRR065654.job;
#	echo "cat SRR065654_fc.bedgraph.tmp | grep -v nan > SRR065654_fc.bedgraph" >> SRR065654.job;
#	qsub SRR065654.job;
#	cp header.job SRR065655.job;
#	echo "macs2 bdgcmp -t SRR065655.bedgraph -c SRR065664.bedgraph -o SRR065655_fc.bedgraph.tmp -m subtract" >> SRR065655.job;
#	echo "cat SRR065655_fc.bedgraph.tmp | grep -v nan > SRR065655_fc.bedgraph" >> SRR065655.job;
#	qsub SRR065655.job;
#	cp header.job SRR065656.job;
#	echo "macs2 bdgcmp -t SRR065656.bedgraph -c SRR065660.bedgraph -o SRR065656_fc.bedgraph.tmp -m subtract" >> SRR065656.job; 
#	echo "cat SRR065656_fc.bedgraph.tmp | grep -v nan > SRR065656_fc.bedgraph" >> SRR065656.job;
#	qsub SRR065656.job;
#	cp header.job SRR065657.job;
#	echo "macs2 bdgcmp -t SRR065657.bedgraph -c SRR065661.bedgraph -o SRR065657_fc.bedgraph.tmp -m subtract" >> SRR065657.job;
#	echo "cat SRR065657_fc.bedgraph.tmp | grep -v nan > SRR065657_fc.bedgraph" >> SRR065657.job;
#	qsub SRR065657.job;
#	cp header.job SRR065658.job;
#	echo "macs2 bdgcmp -t SRR065658.bedgraph -c SRR065662_3.bedgraph -o SRR065658_fc.bedgraph.tmp -m subtract" >> SRR065658.job;
#	echo "cat SRR065658_fc.bedgraph.tmp | grep -v nan > SRR065658_fc.bedgraph" >> SRR065658.job;
#	qsub SRR065658.job;
#	cp header.job SRR065659.job;
#	echo "macs2 bdgcmp -t SRR065659.bedgraph -c SRR065664.bedgraph -o SRR065659_fc.bedgraph.tmp -m subtract" >> SRR065659.job;
#	echo "cat SRR065659_fc.bedgraph.tmp | grep -v nan > SRR065659_fc.bedgraph" >> SRR065659.job;
#	qsub SRR065659.job;

#	   echo "genomeCoverageBed -ibam $(a)_sorted.bam > $(a).bed" >> $(a).job; \
#	echo "genomeCoverageBed -bg -ibam $(a)_filtered.bam > $(a).bedgraph" >> $(a).job; \
#	   echo "samtools mpileup $(a)_sorted.bam -o 

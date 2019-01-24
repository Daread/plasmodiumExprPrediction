include $(HOME)/Templates/Make/quick.mak

CHR = $(THISDIR)
MOTIFSET = $(PARENTDIR)

MOTIFFILE = /net/noble/vol2/home/katecook/proj/2016predictExpression/data/motifs/$(PARENTDIR)/Plasmodium_falciparum.meme
SEQFILE = /net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/chromosomes_Pf/$(CHR).fa

BG = ../../fimo.background

PVALUE_THRESH = 0.01

targets = fimo_out.txt

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)

fimo_out.txt:
	echo "#!/bin/bash" > $@.job;
	echo "#$$ -cwd" >> $@.job;
	echo "#$$ -V" >> $@.job;
	echo "#$$ -l h_vmem=8G" >> $@.job;
	echo "#$$ -l h_rt=36000" >> $@.job;
	echo "source ~/myenv.sh" >> $@.job;
	echo hostname >> $@.job;
	echo "fimo --bgfile $(BG) --text --o $@ --thresh $(PVALUE_THRESH) $(MOTIFFILE) $(SEQFILE) > $@" >> $@.job;
	qsub $@.job;

merge_hits: fimo_out.txt
	echo "#!/bin/bash" > $@.job;
	echo "#$$ -cwd" >> $@.job;
	echo "#$$ -V" >> $@.job;
	echo "#$$ -l h_vmem=8G" >> $@.job;
	echo "#$$ -l h_rt=36000" >> $@.job;
	echo "source ~/myenv.sh" >> $@.job;
	echo hostname >> $@.job;
	echo "python ../../Lib/fimo_to_bedgraph.py $< merged" >> $@.job;
	qsub $@.job;

include $(HOME)/Templates/Make/quick.mak

GENOMEDIR=/net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/chromosomes_Pf
CHRSIZES=/net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/chromosomeLengths

ID = $(THISDIR)

WINDOW_SIZE = 100

CHRLENGTH = $(shell cat $(CHRSIZES) | grep $(ID) | cut -f 2)

targets =

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)


GC.txt: $(GENOMEDIR)/$(ID).fa
	python ../Lib/gc.py $< $@ $(WINDOW_SIZE);

GC.bedgraph: GC.txt
	seq 1 $(CHRLENGTH) \
	| awk 'BEGIN {OFS="\t"}; {print $$1-1,$$1}' \
	| paste - GC.txt \
	| paste.pl $(ID) - \
	> $@;	
	

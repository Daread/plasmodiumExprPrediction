include $(HOME)/Templates/Make/quick.mak

CHILDREN = $(shell cat ids.lst)

GENOME = /net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/PlasmoDB-29_Pfalciparum3D7_Genome.fasta
#GENOME = /net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/chromosomes/chr1.fa

targets = 

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)

echo:
	echo $(CHILDREN)

maker:
	$(foreach c, $(CHILDREN), \
	   mkdir -p $(c); \
	   cd $(c); \
	   ln -sf ../Lib/child.dataset.mak Makefile; \
	   make maker; \
	   cd ..; \
	)


doit:
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   echo $(c); \
	   make assemble_bedgraphs; \
	   cd ..; \
	)

fimo.background:
	fasta-get-markov -m 1 $(GENOME) $@

info.tab:
	$(foreach c, $(CHILDREN), \
	   ls -1 $(c)/*.bedgraph \
	   | cut -f 2 -d '/' \
	   | sed 's/.bedgraph//' \
	   | add_column.pl - -b -s '$(c)' \
	   >> $@; \
	)




include $(HOME)/Templates/Make/quick.mak

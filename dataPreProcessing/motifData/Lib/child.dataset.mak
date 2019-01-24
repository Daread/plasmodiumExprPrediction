include $(HOME)/Templates/Make/quick.mak

ID = $(THISDIR)

SAMPLE_INFO_DIR = /net/noble/vol2/home/katecook/proj/2016predictExpression/results/katecook/sample_info
CHILDREN = $(shell cat $(SAMPLE_INFO_DIR)/chrs_Pf3D7.lst | grep -v PFC10_API_IRAB | grep -v M76611)

TFIDS = $(shell ls */merged/* |  cut -f 3 -d '/' | sed 's/.bedgraph//'  | sort -u)

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
	   ln -sf ../../Lib/child.chr.mak Makefile; \
	   cd ..; \
	)


doit:
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   echo $(c); \
	   make merge_hits; \
	   cd ..; \
	)

assemble_bedgraphs:
	$(foreach t, $(TFIDS), \
	   rm -f $(ID)_$(t).bedgraph; \
	   $(foreach c, $(CHILDREN), \
	      cat $(c)/merged/$(t).bedgraph >> $(ID)_$(t).bedgraph; \
	   ) \
	) 




include $(HOME)/Templates/Make/quick.mak

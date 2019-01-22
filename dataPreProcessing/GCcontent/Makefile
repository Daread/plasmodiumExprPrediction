include $(HOME)/Templates/Make/quick.mak


ID = $(THISDIR)

CHILDREN = $(shell cat ids.lst)

targets =

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)



maker:
	$(foreach c, $(CHILDREN), \
	   mkdir -p $(c); \
	   cd $(c); \
	   ln -sf ../Lib/child.chr.mak Makefile; \
	   make maker; \
	   cd ..; \
	)


doit:
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   echo $(c); \
	   make GC.bedgraph; \
	   cd ..; \
	)

all_GC.bedgraph:
	$(foreach c, $(CHILDREN), \
	   cat $(c)/GC.bedgraph \
	   >> $@; \
	)


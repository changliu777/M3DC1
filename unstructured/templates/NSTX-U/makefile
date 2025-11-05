dirs = fw1_mesh rw1_mesh

.PHONY : $(dirs) clean all
all : $(dirs)

$(dirs) :
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~

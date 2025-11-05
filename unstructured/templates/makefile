dirs := $(wildcard */.)

.PHONY : $(dirs) clean all

all : $(dirs)

$(dirs) : 
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean : $(dirs)
	rm -f *~

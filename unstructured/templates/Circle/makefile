SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)

clean: $(SUBDIRS)
	rm *~

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean $(SUBDIRS)

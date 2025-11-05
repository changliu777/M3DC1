SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)

clean: $(SUBDIRS)
	rm -f *~

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean $(SUBDIRS)

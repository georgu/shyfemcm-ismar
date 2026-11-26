
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FEMDIR = ..
BINDIR = $(FEMDIR)/bin
SUBDIRS = external utils contrib shyfem tools 

MAKEFLAGS += --no-print-directory

#------------------------------------------------------------------------

default: shyfem

fem: shyfem
shyfem: libmod nothing
	@$(BINDIR)/recursivemake fem $(SUBDIRS)

test: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

regress: nothing

status: nothing

list:
	pwd; ls

libmod: $(FEMDIR)/lib/mod
$(FEMDIR)/lib/mod:
	mkdir -p $(FEMDIR)/lib/mod

depend: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

clean: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)
	-rm -f *.tmp *.bak

cleanall: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)
	-rm -f *.tmp *.bak

save: nothing

zip: nothing

nothing:
	@true

#------------------------------------------------------------------------


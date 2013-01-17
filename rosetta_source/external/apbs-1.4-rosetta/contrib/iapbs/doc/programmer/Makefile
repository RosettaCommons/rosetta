#
# $Id: Makefile,v 1.2 2005/05/24 04:48:53 rok Exp $
#

include ../../etc/iapbs-version.mk
SED	= sed
PKGROOT	= /export/html/iapbs/programmer
DIST.ROOT = ./iapbs-programmer-$(VERSION)

SEDSPEC += \
        -e 's%@VERSION@%$(VERSION)%g' \
        -e 's%@PUBDATE@%$(PUBDATE)%g'

html: 
	$(SED) $(SEDSPEC) Doxyfile.in > Doxyfile
	doxygen Doxyfile

pdf: html
	(cd latex; make)

build: html pdf

install:: build
	install -ma+r html/* $(PKGROOT)/
	install -ma+r latex/refman.pdf $(PKGROOT)/iapbs-programmer.pdf

dist:: build
	mkdir -p $(DIST.ROOT)
	install -ma+r html/* $(DIST.ROOT)/
	install -ma+r latex/refman.pdf $(DIST.ROOT)/iapbs-programmer.pdf
	tar cvzf $(DIST.ROOT).tar.gz $(DIST.ROOT)/

clean: 
	rm -rf html/ latex/ Doxyfile $(DIST.ROOT).tar.gz $(DIST.ROOT)



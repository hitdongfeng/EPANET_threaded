#---------------------------------------------------------------------
#
# Makefile of EPANET2 for Linux
#
VERSION = threaded.r177svn
#
# $LastChangedBy: manu $ $Revision: 175 $
# $LastChangedDate: 2008-09-10 16:19:29 +0200 (Wed, 10 Sep 2008) $
#
#---------------------------------------------------------------------
#
# Copyright (c) 2005, 2006 Manuel Lopez-Ibanez
# TeX: \copyright 2005, 2006 Manuel L\'opez-Ib\'a\~nez
#
# This program is free software (software libre); you can redistribute
# it and/or modify it under the terms of version 2 of the GNU General
# Public License version as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, you can obtain a copy of the GNU
# General Public License at:
#                http://www.gnu.org/copyleft/gpl.html
# or by writing to:
#          Free Software Foundation, Inc., 59 Temple Place,
#                Suite 330, Boston, MA 02111-1307 USA
#
#---------------------------------------------------------------------
PWD	:= $(shell /bin/pwd)

DEBUG=1
# =2 extra information is generated and checks are perfomed
# =1 only checks are performed
# =0 disable (fastest).

LIBTOOLKIT=libtoolkit.a
LIB_DIR =../

SRC_DIR=./

OBJS = epanet.o hash.o hydraul.o inpfile.o input1.o input2.o input3.o	\
	mempool.o output.o quality.o report.o rules.o smatrix.o

HEADERS = enumstxt.h funcs.h hash.h mempool.h text.h toolkit.h types.h	\
	 vars.h

OBJS :=$(foreach OBJ, $(OBJS), $(join $(SRC_DIR), $(OBJ)) )
HEADERS :=$(foreach OBJ, $(HEADERS), $(join $(SRC_DIR), $(OBJ)) )


CC = gcc
OPTIMISE := -O0
override OPTIMISE += -finline-functions -Winline -fmerge-constants
ifdef march
override OPTIMISE += -march=$(march)
endif

ifneq ($(DEBUG),0)
CDEBUG = -g
endif
override CFLAGS += $(OPTIMISE) -Wall -W -D DEBUG=$(DEBUG) $(CDEBUG)

DELETE = @rm -f
ECHO = @echo "$(1)"

SVN_REV := $(if $(shell which svnversion &> /dev/null && [ `svnversion -n .` != exported ] && echo 1),$(shell svnversion -n . | tee svn_version),$(shell cat svn_version))

update_version = @if ! grep -q -x -F "define EPANET_VERSION $(VERSION)" $(1);\
	then \
	sed -r -i 's/^(\#define[ ]+EPANET_VERSION).*/\1 \"$(VERSION)\"/' $(1);\
	sed -r -i 's/^(VERSION:) see Makefile/\1 $(VERSION)/' $(1);\
	fi

DIST_SRC := EPANET-$(VERSION)-src
DIST_BIN := EPANET-$(VERSION)-lib

.PHONY : all cleanall cleanexe clean toolkit dist tags

# Create static library.
toolkit : CFLAGS += -D LINUX
toolkit : $(LIB_DIR)/$(LIBTOOLKIT)

$(LIB_DIR)/$(LIBTOOLKIT) : $(OBJS) $(HEADERS)
	$(call ECHO,--> Creating static library \"$(LIB_DIR)/$(LIBTOOLKIT)\" version $(VERSION) <---)
	ar rcs $(LIB_DIR)/$(LIBTOOLKIT) $(OBJS)
	cp toolkit.h $(LIB_DIR)/
	$(call update_version,$(LIB_DIR)/toolkit.h)

$(OBJS) : $(HEADERS)

%.o : %.c
	$(CC)  $(CFLAGS) -c -o $@ $<

cleanall : clean cleanexe

cleanexe :
	$(call ECHO,---> Removing $(LIB_DIR)/$(LIBTOOLKIT) <---)
	$(DELETE) $(LIB_DIR)/$(LIBTOOLKIT)

clean :
	$(call ECHO,---> Removing object files <---)
	$(DELETE) $(OBJS)
	$(DELETE) *~


all :  cleanall toolkit

dist : DEBUG=0
dist : CDEBUG=
dist : all
	@(mkdir -p ../$(DIST_SRC) && rsync -aC . ../$(DIST_SRC)/ && cd .. \
	&& sed -i 's/^\(VERSION *[=:] \).*/\1$(VERSION)/' \
	   $(DIST_SRC)/Makefile $(DIST_SRC)/toolkit.h \
	&& tar cf - $(DIST_SRC) | gzip -f9 > $(DIST_SRC).tar.gz \
	&& rm -f ./$(DIST_SRC)/* && rmdir ./$(DIST_SRC)/ \
	&& echo "$(DIST_SRC).tar.gz created." && cd $(PWD) )
	@(mkdir -p ../$(DIST_BIN) && cd .. \
	&& cp -f toolkit.h libtoolkit.a ./$(DIST_BIN)/ \
	&& tar cf - $(DIST_BIN) | gzip -f9 > $(DIST_BIN).tar.gz \
	&& rm -f ./$(DIST_BIN)/* && rmdir ./$(DIST_BIN)/ \
	&& echo "$(DIST_BIN).tar.gz created." && cd $(PWD))

tags: TAGS
TAGS: $(OBJS:.o=.c) $(HEADERS)
	etags $(OBJS:.o=.c) $(HEADERS)

-include local.mk


####
# inventory
####

BINARIES = lens2mat mat2lens ptmatch reg2pts

ifdef HAVE_CFITSIO
BINARIES += immap srcim
ifdef HAVE_REGIONS
BINARIES += regcrop
endif
endif

ifndef NO_PYTHON
BINARIES += dist
endif


####
# compiler and linker settings
####

# general settings
CFLAGS = -std=c99 -Wall -Werror -pedantic
LDFLAGS = 
LDLIBS = -lm

# debug or release settings
ifdef DEBUG
CFLAGS += -O0 -g -DDEBUG
else
CFLAGS += -Ofast
endif


####
# build rules
####

.PHONY: all clean

all: $(BINARIES)

clean:
	$(RM) $(BINARIES)

immap: src/immap.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lcfitsio

lens2mat: src/lens2mat.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

mat2lens: src/mat2lens.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

ptmatch: src/ptmatch.c src/input.c src/mpfit.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

regcrop: src/regcrop.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lregions -lcfitsio

reg2pts: src/reg2pts.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

srcim: src/srcim.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lcfitsio

%: src/%.py
	cp $< $@ && chmod +x $@

-include local.mk


####
# inventory
####

ALL = lens2mat mat2lens ptmatch reg2pts

ifdef HAVE_CFITSIO
ALL += immap ptcrop srcim
endif

ifndef NO_PYTHON
ALL += dist
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

all: $(ALL)

clean:
	$(RM) $(ALL)

immap: src/immap.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lcfitsio

lens2mat: src/lens2mat.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

mat2lens: src/mat2lens.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

ptcrop: src/ptcrop.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lcfitsio

ptmatch: src/ptmatch.c src/input.c src/mpfit.c src/mpis.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

reg2pts: src/reg2pts.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

srcim: src/srcim.c src/input.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lcfitsio

%: src/%.py
	cp $< $@ && chmod +x $@

####
# inventory
####

BINARIES = ptmap regcrop regpts


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

ptmap: src/ptmap.c src/input.c src/newuoa.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

regcrop: src/regcrop.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS) -lregions -lcfitsio

regpts: src/regpts.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

CC=pgcc

CFLAGS=-lm -m64

SOURCES=mcmc.c em.c

OBJECTS=$(SOURCES:.c=.o)

all: $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@



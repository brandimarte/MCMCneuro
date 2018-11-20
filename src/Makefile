# ***************************************************************
#
# Makefile for gcc compiler
#
# ***************************************************************

CFLAGS = -g -I. -O2 -Wall -pedantic -ansi
LDFLAGS = -O2 -L.
LDLIBS = -lm 

RM = /bin/rm -f
CC = gcc

#======================================================================

.c.o:
	$(CC) $(CFLAGS) -c $*.c

.c:
	make $*.o
	$(CC) $(CFLAGS) -o $* $*.o $(LDLIBS) 

#======================================================================

clean:
	$(RM) *~ \#~ .\#* *.o core a.out

#======================================================================

graphPenalty: Utils.o Item.o ST.o Neuro.o graphPenalty.o
	$(CC) $(CFLAGS) -o ../bin/graphPenalty Utils.o Item.o ST.o Neuro.o graphPenalty.o $(LDLIBS) 

bestGraph: Utils.o Item.o ST.o Neuro.o bestGraph.o
	$(CC) $(CFLAGS) -o ../bin/bestGraph Utils.o Item.o ST.o Neuro.o bestGraph.o $(LDLIBS) 


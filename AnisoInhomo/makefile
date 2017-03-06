#	Makefile for rfgtest
#	Author: A.Smirnov, www.cemr.wvu.edu/~andrei

FC=g77 -c
CC=cc -c
LD=g77

PROG=rfgtest

OBJ=rfgtest.o rfg.o gauss.o random.o diag.o tred2.o tql2.o phytag.o

all:	$(PROG)

$(PROG): $(OBJ)
	$(LD) $(OBJ) -o $(PROG)

rfg.o: rfg.c rfg.h vecalg.h
	$(CC) rfg.c

.c.o:
	$(CC) $<

.f.o:
	$(FC) $<

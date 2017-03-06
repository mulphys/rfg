#	Makefile for rfgtest

FC=g77 -c
CC=cc -c
#LD=g77 # for fortran
LD=cc -lm

PROG=rfgtest

#OBJ=rfgtestf.o rfg.o gauss.o random.o # for fortran
OBJ=rfgtest.o rfg.o gauss.o 

all:	$(PROG)

$(PROG): $(OBJ)
	$(LD) $(OBJ) -o $(PROG)

rfg.o: rfg.c rfg.h vecalg.h
	$(CC) rfg.c

.c.o:
	$(CC) $<

.f.o:
	$(FC) $<

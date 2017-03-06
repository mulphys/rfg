#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vecalg.h"
#include "rfg.h"

#ifndef SMALL
#define SMALL	1.e-30
#endif

int	ne = 1; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

REAL
	*Omega, /* Eq.15,celik.bib:\cite{LiAhetalJAS94} */
	*U1,*U2, /* velocity vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
	*K; /* wave vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */

void    Allocate(REAL **A, int n)
{
	if ((*A = (REAL *) malloc(sizeof(**A)*n)) == NULL)
	{
		fprintf(stderr,"CAN'T ALLOCATE MEMORY\n");
		exit(1);
	}
}
void	genspec_
(
	int	*Ne,
	REAL	*TURB_TIME,
	REAL	*TURB_LENGTH
)
/*
	Generate spectral expansion coefficients 
*/
{
	extern	REAL	gauss_(); /* get a Gaussian random variable */
	extern	void	gaussn_(REAL *, REAL, int); /* get an array of Gaussian random numbers */
	int	ie;
	REAL	fe,turb_time=*TURB_TIME;
	ne=*Ne;
	if (ne<=0) return;
	fe=sqrt(2./(REAL)ne);
	Allocate(&Omega,ne);
	Allocate(&K ,ne*DIM);
	Allocate(&U1,ne*DIM);
	Allocate(&U2,ne*DIM);
	seed_(); /* initialize the random generator */
	gaussn_(K,.5,ne*DIM);
	for (ie=0; ie<ne; ie++)
	{	int	i,j=ie*DIM;
		REAL
			a,V1[DIM],V2[DIM], /* random vectors xi and zeta (Eq.16) */
			*u1=U1+j,*u2=U2+j,*k=K+j;
		Omega[ie] = gauss_()/turb_time; 
		gaussn_(V1, fe, DIM);
		gaussn_(V2, fe, DIM);
		VECP(u1,V1,k); /* Eq.16\cite{LiAhetalJAS94} */
		VECP(u2,V2,k);
		for (i=0; i<DIM; i++)k[i]/=TURB_LENGTH[i];
	}
}
void	delspec_()
{
	free(Omega);
	free(K);
	free(U1);
	free(U2);
}
void	genvec_
(
//	INPUT:
	REAL	*t,  // time
	REAL	*x,  // coordinates
// OUTPUT:
	REAL	*v   // velocities
)
{
	int	i,ie;
	REAL	a,c,s;
	for (i=0; i<DIM; i++) v[i]=0.0;
	if (ne<=0)  return;
	for (ie=0; ie<ne; ie++)
	{	int n=ie*DIM;
		REAL	*k=K+n;
		a=SCLP(k,x)+Omega[ie]**t;
		c=cos(a); s=sin(a);
		for (i=0; i<DIM; i++) 
			v[i]+=U1[n+i]*c+U2[n+i]*s;
	}
}

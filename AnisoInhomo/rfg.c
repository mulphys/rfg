#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rfg.h"
#include "vecalg.h"

#ifndef SMALL
#define SMALL	1.e-30
#endif

int	ne = 1; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

REAL
	fe = 1.41421,
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
	int	*Ne
)
/*
 *	Generate spectral expansion coefficients 
 */
{	extern	void	diag
	(
		REAL *A, /* velocity correlations */  
		REAL *D  /* diagonal vector after diagonalization of A */
	);
	extern	REAL	gauss_(); /* get a Gaussian random variable */
	extern	void	gaussn_(REAL *, REAL, int); /* get an array of Gaussian random numbers */
	int	i,ie;
	ne=*Ne;
	if (ne<=0) return;
	fe=sqrt(2./(REAL)ne);
	Allocate(&Omega,ne);
	Allocate(&K ,ne*DIM);
	Allocate(&U1,ne*DIM);
	Allocate(&U2,ne*DIM);
	seed_(); /* initialize the random generator */
//	seed0(999); //same numbers every time
	gaussn_(K,.5,ne*DIM);
	for (ie=0; ie<ne; ie++)
	{	int	j=ie*DIM;
		REAL	a,V1[DIM],V2[DIM], /* random vectors xi and zeta (Eq.16) */
			*u1=U1+j,
			*u2=U2+j,
			*k=K+j;
		Omega[ie] = gauss_(); 
		gaussn_(V1, 1., DIM);
		gaussn_(V2, 1., DIM);
		VECP(u1,V1,k); /* Eq.16\cite{LiAhetalJAS94} */
		VECP(u2,V2,k);
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
	REAL	*TT, // Turbulent Time: scalar
	REAL	*UU, // velocity correlations: UU,UV,VV,UW,VW,WW
// OUTPUT:
	REAL	*v   // velocities
)
{	extern	void	diag
	(
		REAL *A, /* velocity correlations */  
		REAL *D  /* diagonal vector after diagonalization of A */
	);
	int	i,ie;
	REAL	a,c,s,
		d[DIM],dd=0.0,
		turb_time=*TT;

	if (ne<=0)  return;
	diag(UU,d);
	for (i=0; i<DIM; i++)
	{	REAL	r=fabs(d[i]);
		d[i]=sqrt(r);
		dd+=r;
		v[i]=0.0;
	}
	dd=sqrt(dd);
	for (ie=0; ie<ne; ie++)
	{	int n=ie*DIM;
		REAL	k[DIM],
		      *u1=U1+n,*u2=U2+n;
		for (i=0; i<DIM; i++)
			k[i]=K[n+i]/(d[i]>SMALL?turb_time*d[i]:turb_time*dd);
		a=SCLP(k,x)+Omega[ie]**t/turb_time;
		c=cos(a); s=sin(a);
		for (i=0; i<DIM; i++) 
			v[i]+=u1[i]*c+u2[i]*s;
	}
	for (i=0; i<DIM; i++)v[i]*=fe*d[i];//V-anisotropy
	btrans(v);
}


/*
	Vector algebra
*/

#ifndef	REAL
#define	REAL	double
#endif
#ifndef	DIM
#define	DIM	3
#endif

extern int	E[][][]; /* assymmetric tensor used in vector product */

#define MUL0(x,s)	for (i=0; i<DIM; i++) (x)[i] *= (s)
#define ADD0(x,y)	for (i=0; i<DIM; i++) (x)[i] += (y)[i]
#define ZERO(x)	for (i=0; i<DIM; i++) (x)[i] = 0.
#define SCLP(A,B)	(A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define VECP(A,B,C)	A[0]=B[1]*C[2]-B[2]*C[1];\
							A[1]=B[2]*C[0]-B[0]*C[2];\
							A[2]=B[0]*C[1]-B[1]*C[0];
#define SCLP1(A,B)	*(A+0)*B[0]+*(A+1)*B[1]+*(A+2)*B[2]
#define LENGTH(A)	sqrt(SCLP(A,A))

void	inivec();
void	vecp
/*
	Vector product A = [B,C]
*/
(
	REAL *A,
	REAL *B,
	REAL *C
);
REAL	average
(
	REAL	*A
);
REAL	sclp
/* 	Scalar product: a = (B,C)
*/
(
	REAL *A,
	REAL *B
);

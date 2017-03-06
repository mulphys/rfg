#include <stdlib.h>
#include <stdio.h>
#include "vecalg.h"

extern void
	genspec_(int*,REAL*,REAL*),
	genvel_(REAL*, REAL*, REAL*),
	delspec_();

	int nspec=1000;     // spectral sample size
double 
	turbtime=1.0,   // turbulence time-scale
	turblength[] = {1.0, 1.0, 1.0};

int main(int argc, char *argv[]) {
	int n; //number of time steps
	REAL
		t = 0.0, dt = 1.0,
		V00, V01, V02, // correlations
		X[DIM],
		V[DIM];

	if (argc != 2) {
		fprintf(stderr,"Usage:\n\t%s <number_of_time_steps>\n",argv[0]);
		return 1;
	}
	n = atoi(argv[1]);
	genspec_(&nspec,&turbtime,turblength);
	X[0] = 1.0;
	X[1] = 1.0;
	X[2] = 1.0;
	V00 = V01 = V02 = 0.0; 
	for (int i=0; i<n; i++,t+=dt) {
		genvel_(&t,X,V);
		//printf("%d: %g\t%g\t%g\n",i+1,V[0],V[1],V[2]);
		V00 += V[0]*V[0];
		V01 += V[0]*V[1];
		V02 += V[0]*V[2];
	}
	V00 /= n;
	V01 /= n;
	V02 /= n;
	printf("Stress tensor averaged over %g realizations:\n",n);
	printf("<Vx*Vx>: %g\n",V00);
	printf("<Vx*Vy>: %g\n",V01);
	printf("<Vx*Vz>: %g\n",V02);
	delspec_();
	return 0;
}

#ifndef __simplex_noise_h__
#define __simplex_noise_h__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>

class simplex {
public:
	simplex(int, const char * = "");
	~simplex();
	double noise(double *);
	double octavednoise(double *, unsigned, double *, double);
	void save(const char *);

private:
	int dim, numgrads;
	int * perm, * iijj, * ranked, * skewed, * gi;
	int ** grads, ** ij;

	double skewing, unskewing, scalar;
	double * dist, * oct;
	double ** xy;

	void calcgrads(int **, int);
	void calcperms(int **, int, int &, int, int, int);

	template<class t1> double sum(t1 *, unsigned);
	template<class t1, class t2> double dot(t1 *, t2 *, unsigned);
};

#endif

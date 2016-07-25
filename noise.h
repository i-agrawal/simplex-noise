#ifndef __noise_h__
#define __noise_h__

#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>

class noise {
public:
	noise(int);
	noise(int, const char *);
	double simplexnoise(double *);
	void save(const char *);

private:
	int dim, numgrads;
	double skewing, unskewing;
	int * perm;
	int ** grads;

	void calcgrads(int **, int);
	void calcperms(int **, int, int &, int, int, int);

	int * rank(double *, unsigned);
	template<class t1> double sum(t1 *, unsigned);
	template<class t1> int * floor(t1 *, unsigned);
	template<class t1> double * constant(t1, unsigned);
	template<class t1, class t2> double dot(t1 *, t2 *, unsigned);
	template<class t1, class t2> double * add(t1 *, t2 *, unsigned);
	template<class t1, class t2> double * add(t1 *, t2, unsigned);
	template<class t1, class t2> double * sub(t1 *, t2 *, unsigned);
	template<class t1, class t2> double * mult(t1 *, t2, unsigned);
	template<class t1> void store(t1 *, const char *);
	template<class t1> void load(t1 *, const char *);
};

#endif

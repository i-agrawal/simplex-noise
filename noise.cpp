#include "noise.h"

noise::noise(int dimensions) {
	//record our number of dimensions
	dim = dimensions;
	int i;

	//begin creating permutation array through rand()
	//initialize permutation array
	perm = (int *)malloc(512 * sizeof(int));

	for (i = 0; i < 256; i++) {
		perm[i] = -1;
	}

	int randindex = 0;
	for (i = 0; i < 256; i++) {
		do { randindex = rand() % 256; }
		while(perm[randindex] != -1);
		perm[randindex] = perm[randindex + 256] = i;
	}
	//end creating permutation array through rand()

	//begin creating 2D gradient array through dimensions
	int perms = pow(2, dim-1);
	int combs = numgrads = dim * perms;
	grads = (int **)malloc(combs * sizeof(int *));
	for (i = 0; i < combs; i++) {
		grads[i] = (int *)malloc(dim * sizeof(int));
	}
	calcgrads(grads, dim);
	//begin creating 2D gradient array through dimensions

	//precalculate skewing and unskewing doubles
	skewing = (sqrt(dim + 1.0) - 1.0) / dim;
	unskewing = (1.0 / (sqrt(dim + 1.0) - 1.0)) / dim;
	//end calculations
}

noise::noise(int dimensions, const char * fn) {
	int i;
	dim = dimensions;

	//begin loading permutation arrary through file
	perm = (int *)malloc(512 * sizeof(int));
	int fd = open(fn, O_RDONLY);
	read(fd, perm, sizeof(perm));
	//end loading of permutation array through file

	//begin creating 2D gradient array through dimensions
	int perms = pow(2, dim-1);
	int combs = numgrads = dim * perms;
	grads = (int **)malloc(combs * sizeof(int *));
	for (i = 0; i < combs; i++) {
		grads[i] = (int *)malloc(dim * sizeof(int));
	}
	calcgrads(grads, dim);
	//end creating 2D gradient array through dimensions
}

void noise::calcgrads(int ** grads, int dim) {
	int i, yindex = 0;
	for (i = 0; i < dim; i++) {
		calcperms(grads, 0, yindex, dim-i-1, i, 1);
	}
}

void noise::calcperms(int ** grads, int x, int &y, int ones, int nones, int zeros) {
	int i;
	int j = ones + nones + zeros;
	bool didother = 0;
	if (ones != 0) {
		grads[y][x] = 1;
		calcperms(grads, x+1, y, ones-1, nones, zeros);
		if (j == 1) {
			y = y+1;
		}
		didother = 1;
	}
	if (nones != 0) {
		grads[y][x] = -1;
		if (didother) {
			for (i = 0; i < x; i++) {
				grads[y][i] = grads[y-1][i];
			}
		}
		calcperms(grads, x+1, y, ones, nones-1, zeros);
		if (j == 1) {
			y = y+1;
		}
		didother = 1;
	}
	if (zeros != 0) {
		grads[y][x] = 0;
		if (didother) {
			for (i = 0; i < x; i++) {
				grads[y][i] = grads[y-1][i];
			}
		}
		calcperms(grads, x+1, y, ones, nones, zeros-1);
		if (j == 1) {
			y = y+1;
		}
	}
}

double noise::simplexnoise(double * point) {
	int i, j, k;
	double shairy = skewing * sum(point, dim);

	int * skewed = floor(add(point, shairy, dim), dim);
	double uhairy = unskewing * sum(skewed, dim);
	double * dist = add(point, add(skewed, -uhairy, dim), dim);

	int * ranked = rank(dist, dim);
	int ** ij = (int **)malloc((dim+1) * sizeof(int *));

	for (i = 0; i < dim+1; i++) {
		ij[i] = (int *)malloc(dim * sizeof(int));
		for (j = 0; j < dim; j++) {
			ij[i][j] = (ranked[j] >= (dim - i + 1));
		}
	}

	double ** xy = (double **)malloc((dim+1) * sizeof(double *));
	xy[0] = dist;
	for (i = 1; i < dim+1; i++) {
		xy[i] = add(sub(dist, ij[i], dim), (i)*unskewing, dim);
	}

	int * iijj = (int *)malloc(dim * sizeof(int));
	for (i = 0; i < dim; i++) {
		iijj[i] = skewed[i] & 255;
	}

	int * gi = (int *)malloc((dim + 1) * sizeof(int));
	for (i = 0; i < dim+1; i++) {
		k = 0;
		for (j = dim-1; j >= 0; j--) {
			k = perm[k + iijj[j] + ij[i][j]];
		}
		gi[i] = k % numgrads;
	}

	double t, n = 0;
	for (i = 0; i < dim+1; i++) {
		t = 0.5 - dot(xy[i], xy[i], dim);
		if (t > 0) {
			n += pow(t, 4) * dot(grads[gi[i]], xy[i], dim);
		}
	}

	double scalar = 1;
	if (dim == 2) {
		scalar = 105.6; //max without = 0.009389
	}
	else if (dim == 3) {

	}

	return scalar * n;
}

int * noise::rank(double * a, unsigned sz) {
	int * ret = (int *)malloc(sz * sizeof(int));
	int i, j, k;
	for (i = 0; i < sz; i++) {
		k = 0;
		for (j = 0; j < sz; j++) {
			if (a[i] > a[j] || (i <= j && a[i] >= a[j])) {
				k++;
			}
		}
		ret[i] = k;
	}
	return ret;
}

template<class t1> double noise::sum(t1 * a, unsigned sz) {
	int i;
	double sum = 0;
	for (i = 0; i < sz; i++) {
		sum += a[i];
	}
	return sum;
}

template<class t1> int * noise::floor(t1 * a, unsigned sz) {
	int * ret = (int *)malloc(sz * sizeof(int));
	int i;
	for (i = 0; i < sz; i++) {
		ret[i] = (a[i] > 0) ? static_cast<int>(a[i]) : static_cast<int>(a[i])-1;
	}
	return ret;
}

template<class t1> double * noise::constant(t1 a, unsigned sz) {
	double * ret = (double *)malloc(sz * sizeof(double));
	int i;
	for (i = 0; i < sz; i++) {
		ret[i] = a;
	}
	return ret;
}

template<class t1, class t2> double noise::dot(t1 * a, t2 * b, unsigned sz) {
	int i;
	double sum = 0;
	for (i = 0; i < sz; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

template<class t1, class t2> double * noise::add(t1 * a, t2 * b, unsigned sz) {
	double * ret = (double *)malloc(sz * sizeof(double));
	int i;
	for (i = 0; i < sz; i++) {
		ret[i] = a[i] + b[i];
	}
	return ret;
}

template<class t1, class t2> double * noise::add(t1 * a, t2 b, unsigned sz) {
	double * ret = (double *)malloc(sz * sizeof(double));
	int i;
	for (i = 0; i < sz; i++) {
		ret[i] = a[i] + b;
	}
	return ret;
}

template<class t1, class t2> double * noise::sub(t1 * a, t2 * b, unsigned sz) {
	double * ret = (double *)malloc(sz * sizeof(double));
	int i;
	for (i = 0; i < sz; i++) {
		ret[i] = a[i] - b[i];
	}
	return ret;
}

template<class t1, class t2> double * noise::mult(t1 * a, t2 b, unsigned sz) {
	double * ret = (double *)malloc(sz * sizeof(double));
	int i;
	for (i = 0; i < sz; i++) {
		ret[i] = a[i] * b;
	}
	return ret;
}

void noise::save(const char * fn) {
	int fd = open(fn, O_WRONLY);
	write(fd, perm, sizeof(perm));
}

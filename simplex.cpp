#include "simplex.h"

simplex::simplex(int dimensions, const char * fn) {
	//begin creation of simplex noise
	int i, randindex;
	dim = dimensions;

	//begin creation permutation array
	perm = (int *)malloc(512 * sizeof(int));
	if (!strcmp(fn, "")) {
		//begin loading of permutation array through file
		int fd = open(fn, O_RDONLY);
		read(fd, perm, sizeof(perm));
		//end loading of permutation array through file
	}
	else {
		//begin creation of permutation array through rand()
		for (i = 0; i < 256; i++) {
			perm[i] = -1;
		}
		for (i = 0; i < 256; i++) {
			do { randindex = rand() % 256; }
			while(perm[randindex] != -1);
			perm[randindex] = perm[randindex + 256] = i;
		}
		//end creation of permutation array through rand()
	}
	//end creation of permutation array

	//begin creation of 2D gradient matrix
	int perms = pow(2, dim-1);
	int combs = numgrads = dim * perms;
	grads = (int **)malloc(combs * sizeof(int *));
	for (i = 0; i < combs; i++) {
		grads[i] = (int *)malloc(dim * sizeof(int));
	}
	calcgrads(grads, dim);
	//end creation 2D gradient matrix

	//begin precalculation of skewing and unskewing doubles
	skewing = (sqrt(dim + 1.0) - 1.0) / dim;
	unskewing = (1.0 / (sqrt(dim + 1.0) - 1.0)) / dim;
	//end precalculation of skewing and unskewing doubles

	//begin creation of arrays to hold noise values
	ij = (int **)malloc((dim + 1)*sizeof(int *));
	xy = (double **)malloc((dim+1) * sizeof(double *));
	for (i = 0; i < dim + 1; i++) {
		ij[i] = (int *)malloc(dim * sizeof(int));
		xy[i] = (double *)malloc(dim * sizeof(double));
	}
	iijj = (int *)malloc(dim * sizeof(int));
	skewed = (int *)malloc(dim * sizeof(int));
	ranked = (int *)malloc(dim * sizeof(int));
	dist = (double *)malloc(dim * sizeof(double));
	oct = (double *)malloc(dim * sizeof(double));
	gi = (int *)malloc((dim + 1) * sizeof(int));
	//end creation of arrays to hold noise values

	//begin calculating the scalar value per dimension
	scalar = 1;
	if (dim == 2) {
		scalar = 105.6; //max without = 0.009389
	}
	else if (dim == 3) {

	}
	//end calculating the scalar value per dimension
	//end creation of simplex noise
}

simplex::~simplex() {
	int i;
	free(perm);
	free(iijj);
	free(gi);
	free(ranked);
	free(skewed);
	free(dist);
	free(oct);
	for (i = 0; i < dim + 1; i++) {
		free(ij[i]);
		free(xy[i]);
	}
	free(ij);
	free(xy);
	for (i = 0; i < numgrads; i++) {
		free(grads[i]);
	}
	free(grads);
}

void simplex::calcgrads(int ** grads, int dim) {
	int i, yindex = 0;
	for (i = 0; i < dim; i++) {
		calcperms(grads, 0, yindex, dim-i-1, i, 1);
	}
}

void simplex::calcperms(int ** grads, int x, int &y, int ones, int nones, int zeros) {
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

double simplex::noise(double * point) {
	//begin simplex noise calculation
	int i, j, k;																//create constants for loops
	double t;																	//create constant for values

	double shairy = skewing * sum(point, dim); 									//calculate skewing hairy factor = skewing factor
																				// 								   * point's sum

	//begin calculating skewed point
	for (i = 0; i < dim; i++) {													//go through all the point values and skew them
		t = point[i] + shairy;													//get the skewed point by adding the skewing hairy factor
		skewed[i] = (t > 0) ? static_cast<int>(t) : static_cast<int>(t)-1; 		//then floor the answer
	}
	//end calculating skewed point

	double uhairy = unskewing * sum(skewed, dim); 								//calculate unskewing hairy factor = unskewing factor
																				// 									* skewed point's sum

	//begin calculating distance
	for (i = 0; i < dim; i++) {													//go through values in the point
		dist[i] = point[i] + skewed[i] - uhairy;								//the distance = point + skewed, floored, then unskewed point
	}
	//end calculating distance

	//begin finding what has the largest distance, second largest distance, etc.
	for (i = 0; i < dim; i++) {													//go through all the distance
		k = 0;																	//set the number of smaller distances to 0
		for (j = 0; j < dim; j++) { 											//go through all the other distances
			if (dist[i] > dist[j] || (i <= j && dist[i] >= dist[j])) {  		//if they are smaller
				k++;															//add 1 to number of smaller distances
			}
		}
		ranked[i] = k;															//record how many distances this distance is larger than
	}
	//end finding what has the largest distance, second largest distance, etc.

	//begin creating dimensional array of
	for (i = 0; i < dim+1; i++) {												//
		for (j = 0; j < dim; j++) {
			ij[i][j] = (ranked[j] >= (dim - i + 1));
		}
	}

	xy[0] = dist;
	for (i = 1; i < dim+1; i++) {
		for (j = 0; j < dim; j++) {
			xy[i][j] = dist[j] - ij[i][j] + i * unskewing;
		}
	}

	for (i = 0; i < dim; i++) {
		iijj[i] = skewed[i] & 255;
	}

	for (i = 0; i < dim+1; i++) {
		k = 0;
		for (j = dim-1; j >= 0; j--) {
			k = perm[k + iijj[j] + ij[i][j]];
		}
		gi[i] = k % numgrads;
	}

	double n = 0;
	for (i = 0; i < dim+1; i++) {
		t = 0.5 - dot(xy[i], xy[i], dim);
		if (t > 0) {
			n += pow(t, 4) * dot(grads[gi[i]], xy[i], dim);
		}
	}

	return (scalar * n + 1.0) / 2.0;
	//end simplex noise calculation
}

template<class t1> double simplex::sum(t1 * a, unsigned sz) {
	int i;
	double sum = 0;
	for (i = 0; i < sz; i++) {
		sum += a[i];
	}
	return sum;
}

template<class t1, class t2> double simplex::dot(t1 * a, t2 * b, unsigned sz) {
	int i;
	double sum = 0;
	for (i = 0; i < sz; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

void simplex::save(const char * fn) {
	int fd = open(fn, O_WRONLY);
	write(fd, perm, sizeof(perm));
}

double simplex::octavednoise(double * point, unsigned octaves, double * freq, double e) {
	if (e == 0)
		return 1;
	int i, j;
	double t;
	double s = 0;
	for (i = 0; i < octaves; i++) {
		t = pow(2, i);
		for (j = 0; j < dim; j++) {
			oct[j] = point[j] + t;
		}
		s += freq[i] * noise(oct);
	}
	s /= sum(freq, octaves);
	s = pow(s, e);
	return s;
}

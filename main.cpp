#include "simplex.h"
#include <time.h>
#include <cmath>
#include <stdio.h>

int main(int argc, char ** argv) {
	int dim = 2;

	double t;
	printf("Creating New 2D Simplex Noise Generator\n");
	t = clock();
	simplex * newgen = new simplex(dim);
	t = clock() - t;
	printf("Time Elapsed | %f sec (clock creation = 1e-6 sec)\n\n", t / CLOCKS_PER_SEC);

	printf("Loading 2D Simplex Noise Generator From File\n");
	t = clock();
	simplex * loadgen = new simplex(dim, "perm.txt");
	t = clock() - t;
	printf("Time Elapsed | %f sec (clock creation = 1e-6 sec)\n\n", t / CLOCKS_PER_SEC);

	double val;
	int i, j;
	double w = 1000, h = 1000;
	double * point = (double *)malloc(dim * sizeof(double));
	point[0] = 30.0 / w - 0.5;
	point[1] = 50.0 / h - 0.5;
	printf("Performing Noise Function On Point 30, 50 in a %g x %g map\n", w, h);
	t = clock();
	val = loadgen->noise(point);
	t = clock() - t;
	printf("Time Elapsed | %f sec (clock creation = 1e-6 sec)\n\n", t / CLOCKS_PER_SEC);


	double ** map = (double **)malloc(h * sizeof(double *));
	for (i = 0; i < h; i++) {
		map[i] = (double *)malloc(w * sizeof(double));
	}
	printf("Performing Simplex Noise Map Generation | (%g x %g) Map\n", w, h);
	t = clock();
	for (i = 0; i < h; i++) {
		for (j = 0; j < w; j++) {
			point[0] = j / w - 0.5;
			point[1] = i / h - 0.5;
			map[i][j] = loadgen->noise(point);
		}
	}
	t = clock() - t;
	printf("Time Elapsed | %f sec (clock creation = 1e-6 sec)\n\n", t / CLOCKS_PER_SEC);

	unsigned octaves = 4;
	double * freq = (double *)malloc(octaves * sizeof(double));
	freq[0] = 1.0;
	freq[1] = 0.5;
	freq[2] = 0.25;
	freq[3] = 0.125;
	double e = 1;

	printf("Performing Octaved Simplex Noise Map Generation | 4 Octaves | (%g x %g) Map\n", w, h);
	t = clock();
	for (i = 0; i < h; i++) {
		for (j = 0; j < w; j++) {
			point[0] = j / w - 0.5;
			point[1] = i / h - 0.5;
			map[i][j] = loadgen->octavednoise(point, octaves, freq, e);
		}
	}
	t = clock() - t;
	printf("Time Elapsed | %f sec (clock creation = 1e-6 sec)\n\n", t / CLOCKS_PER_SEC);

	return 0;
}

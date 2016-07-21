#include <iostream>
#include <Eigen/Dense>
#include "test.cpp"

typedef Eigen::MatrixXd mat;

int * permcreate() {
	int * perm = new int[512];
	int * p = new int[256];
	int * p2 = new int[256];
	int rng;
	for (unsigned i = 0; i < 256; i++) {
		rng = rand() % 256;
		while(p2[rng] == 1) {
			rng = rand() % 256;
		}
		p[i] = rng;
		p2[rng] = 1;
	}

	for (int i = 0; i < 512; i++)
		perm[i] = p[i & 255];

	return perm;
}

double octavednoise(const mat &input, int * perm, int octaves, double expo) {
	int dim = input.cols();
	mat div = mat::Constant(octaves, 1, 1);
	mat mult = mat::Constant(octaves, 1, 1);
	for (unsigned i = 0; i < octaves; i++) {
		div(i, 0) = 1.0 / pow(2, i); 
		mult(i, 0) = 1.0 * pow(2, i);
	}

	mat a = mult * input;
	mat b(octaves, 1);
	for (unsigned i = 0; i < octaves; i++) {
		b(i, 0) = div(i, 0) * noise(a.row(i), perm);
	}
	double ret = b.sum() / div.sum();
	ret = pow(ret, expo);
	return ret;
}

int main() {
	int * a = permcreate();
	int * b = permcreate();

	int h = 10, w = 10;
	mat hmap(h, w);
	mat mmap(h, w);
	mat i(1, 2);

	for (double y = 0; y < h; y++) {
		for(double x = 0; x < h; x++) {
			double nx = x/w - 0.5, ny = y/h - 0.5;
			i << nx, ny;
			hmap(y, x) = octavednoise(i, a, 6, 1);
			mmap(y, x) = octavednoise(i, b, 6, 1);
		}
	}

	std::cout << hmap << std::endl;

	return 0;
}
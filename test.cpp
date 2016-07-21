#include <iostream>
#include <random>
#include <Eigen/Dense>

typedef Eigen::MatrixXd mat;

//creates a matrix with all the combinations of 1, -1, and a single 0
mat combo(int n) {
	int rs = pow(2, n - 1);
    mat a(n * rs, n);

    for (unsigned i = 0; i < n; i++) {
        int nones = n - 1;
        mat indices = mat::Constant(1, n - 1, 1);
        for (unsigned j = 0; j < rs; j++) {
            a(i * rs + j, i) = 0;
            unsigned l = 0;
            unsigned k = 0;
            while (k < n) {
                if(k != i) {
                    a(i * rs + j, k) = (indices(0, l) == 1) ? 1 : -1;
                    l++;
                }
                k++;
            }
            int ones = nones;
            int ind = n - 2;
            while(ind >= 0 && indices(0, ind) == 1) {
                ones--;
                ind--;
            }
            if (ones == 0 && nones != 0) {
                nones--;
                indices = mat::Constant(1, n - 1, 0);
                for (unsigned m = 0; m < nones; m++) {
                    indices(0, m) = 1;
                }
            }
            else if(nones != 0){
                ind = n-2;
                while(ind >= 0 && indices(0, ind) == 1) {
                    ind--;
                }
                while(ind >= 0 && indices(0, ind) == 0) {
                    ind--;
                }
                indices(0, ind + 1) = 1;
                indices(0, ind) = 0;
            }
        }
    }

    return a;
}

//input is the point in space and dim is the dimensions of the calculation
double noise(const mat &input, int perm[]) {

	int dim = input.cols();
	double F = (sqrt(dim + 1) - 1) / dim;			//skewing factor is (sqrt(dimension + 1) - 1) / dimensions
	double S = F * input.sum();						//skewing hairy factor is the skewing factor times the sum of the inputs
	double G = -1 * (1 / (sqrt(dim + 1)) - 1) / dim; 	//the unskew factor is equal to (1 / sqrt(dimension + 1) - 1) / dimension

	mat skewed = input + mat::Constant(1, dim, S);	//get the vector of skewed coordinateds using hairy factor
	skewed = skewed.cast<int>().cast<double>();		//floor the answers
	double T = G * skewed.sum();					//unskewing hairy factor is the sum of the skewed inputs times the unskewing factor

	mat dist = input - (skewed - mat::Constant(1, dim, T)); //get the distance between the inputs and the unskewed points
															//skewed - unskewing hairy factor unskews the points
	mat simplex(dim, dim);							//this will figure out what simplex we are in

	//we need to figure out what has the largest distance so we are going to rank smallest to largest
	mat asc = dist;
	std::sort(asc.data(), asc.data() + asc.size());
	mat ranked(1, dim);
	for (unsigned i = 0; i < asc.cols(); i++) {
		for (unsigned j = 0; j < dist.cols(); j++) {
			if(asc(0, i) == dist(0, j)) {
				ranked(0, i) = j;
			}
		}
	}

	//then we put a 1 on the highest one in the 1st row and a 1 in the highest 2 in the 2nd row and ...
	for (unsigned i = 0; i < dim; i++) {
		mat temp = ranked;
		for (unsigned j = 0; j < dim; j++) {
			temp(0, j) = (temp(0, j) >= (dim - i - 1))? 1:0;
		}
		simplex.row(i) = temp;
	}
	//now we have the simplex matrix which holds what type of simplex model this is

	//get the offset for each corner
	mat offset(dim, dim);
	for (unsigned i = 0; i < dim; i++) {
		offset.row(i) = dist - simplex.row(i) + mat::Constant(1, dim, (i+1) * G);
	}

	//calculating hashed gradient indices first calculated index in the perm array
	mat hashes(1, dim);
	for (unsigned i = 0; i < dim; i++) {
		hashes(0, i) = (int)(skewed(0, i)) & 255;
	}

	mat grad = combo(dim);

	//this is for utility
	mat simplex2(dim + 1, dim);
	simplex2 << mat::Constant(1, dim, 0),
				simplex;

	//calculate the gradient indices
	mat gindices(1, dim + 1);
	int place = 0;
	for (unsigned i = 0; i < dim + 1; i++) {
		place = 0;
		for(int j = dim - 1; j >= 0; j--) {
			place += hashes(0, j) + simplex2(i, j);
			place = perm[(int)place];
		}
		gindices(0, i) = place % grad.rows();
	}

	mat offset2(dim + 1, dim);
	offset2 << 	dist,
				offset;

	mat ncont(1, dim + 1); 									//there are dimensions + 1 noise contributions to be calculated later
	mat sqdiff(1, dim + 1);
	for (unsigned i = 0; i < dim + 1; i++) {
		sqdiff(0, i) = 0.5;
		for (unsigned j = 0; j < dim; j++) {
			sqdiff(0, i) = sqdiff(0, i) - offset2(i, j) * offset2(i, j);
		}
	}

	for (unsigned i = 0; i < dim + 1; i++) {
		if(sqdiff(0, i) <= 0)
			ncont(0, i) = 0;
		else {
			ncont(0, i) = sqdiff(0, i) * sqdiff(0, i) * sqdiff(0, i) * sqdiff(0, i) * grad.row(gindices(0, i)).dot(offset2.row(i));
		}
	}

	int mscalar = 0;
	if(dim == 2) {
		mscalar = 70;
	}
	else if(dim == 3) {
		mscalar = 32;
	}
	else {
		mscalar = 27;
	}

	double ret = mscalar * ncont.sum();
	return ret / 2.0 + 0.5;

	//return mscalar * ncont.sum();
}
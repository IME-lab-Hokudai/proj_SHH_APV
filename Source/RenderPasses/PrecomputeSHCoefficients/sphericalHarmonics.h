#pragma once
#include <vector>
#include <math.h>
#include <iostream>
using namespace std;

#ifndef M_PI
#define M_PI 3.14159268
#endif

void tchebycheff(
	vector<double>& p,
	const double x,
	const int l
);
void legendre(
	vector<double>& p,
	const double x,
	const int l
);
void alegendre(
	vector<double>& p,
	const double x,
	const int l,
	const int m
);
void calcNormalizationCoeff(
	vector<double>& coeff,
	const int degree
);
void sphericalHarmonics(
	vector<double>& sh,
	const vector<double>& normalizedCoeffs,
	const int degree,
	const double ct,
	const double phi
);
void sphericalHarmonics(
	vector<double>& sh,
	const int degree,
	const double ct,
	const double phi
);

//============================================
class SphericalHarmonics {
protected:
	int		degree;
	int		numBasis;
	vector<double> coeff;
public:

	SphericalHarmonics(int _degree = 8) {
		degree = _degree;
		numBasis = (_degree+1)*(_degree+1);
		calcNormalizationCoeff(coeff, degree);
	}
	~SphericalHarmonics() {}

	int getDegree() {
		return(degree);
	}
	int getNumBasis() {
		return(numBasis);
	}
        //calculate final numerical value of an SH basis
       //ct is cos(theta)
	void calcSHBasis(
		vector<double>& out,
		double ct,
		double phi
	) {
            if (out.size() != numBasis)
                out.resize(numBasis);
            sphericalHarmonics(out, coeff, degree, ct, phi);
	}
};

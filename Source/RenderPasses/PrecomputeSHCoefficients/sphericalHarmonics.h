#pragma once
#include <vector>
#include <math.h>
#include <iostream>

#include "Utils/Math/VectorTypes.h"
using namespace std;
using namespace Falcor;
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
	vector<double>& shBasis,
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

// Helper: Condon-Shortley Phase Factor (-1)^m
// Since std::assoc_legendre OMITS this phase, we must multiply by it
// if the target basis requires it (which your mismatch suggests).
double CSPhase(int m);

/**
 * Computes the normalization constant K_l^m for Spherical Harmonics.
 * Based on the standard SH definition implied by the paper's references.
 */
double getNormalizationK(int l, int m);

/**
 * Computes the derivative of the Associated Legendre Polynomial dP/dx
 * using the recurrence formula from Eq 708 in the paper.
 * * Formula: dP/dx = (1 / (x^2 - 1)) * (x * P_l^m - (l + m) * P_{l-1}^m)
 */
double get_dP_dx(int l, int m, double x);

/**
 * Computes the partial derivative of the Spherical Harmonic y_l^m
 * with respect to theta.
 * * Implements the piecewise formula from Appendix B Krivanek 2005
 *
 * @param l      Band index (degree)
 * @param m      Mode index (order)
 * @param theta  Polar angle in radians [0, pi]
 * @param phi    Azimuthal angle in radians [0, 2pi)
 */
double computeThetaGradient(int l, int m, double theta, double phi);

/**
 * Computes the partial derivative of the Spherical Harmonic y_l^m
 * with respect to phi.
 * * Formula from Appendix B:
 * dy/dphi = -m * y_l^{-m}(theta, phi)
 */
double computePhiGradient(int l, int m, double ylmminus);


//============================================
class SphericalHarmonics {
protected:
	int		degree;
	int		numBasis;
	vector<double> normalizedCoeffs;
public:

	SphericalHarmonics(int _degree = 8) {
		degree = _degree;
		numBasis = (_degree+1)*(_degree+1);
		calcNormalizationCoeff(normalizedCoeffs, degree);
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
            sphericalHarmonics(out, normalizedCoeffs, degree, ct, phi);
	}
};

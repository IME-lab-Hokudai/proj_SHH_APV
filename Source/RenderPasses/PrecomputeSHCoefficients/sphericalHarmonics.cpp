#include "sphericalHarmonics.h"

//-----------------------------------------------------------------------
void tchebycheff(
	vector<double>& p,
	const double x,
	const int l
)
{
	p.resize(l + 1);
	p[0] = 1.0;
	p[1] = x;
	for (int i = 2; i <= l; i++) {
		p[i] = 2.0 * x * p[i - 1] - p[i - 2];
	}
}
//-----------------------------------------------------------------------
void legendre(
	vector<double>& p,
	const double x,
	const int l
)
{
	p.resize(l + 1);
	p[0] = 1.0;
	p[1] = x;
	for(int i=2; i<=l; i++) {
		p[i] = (2.0*(double)i-1.0)*x*p[i-1]-(double)(i-1.0)*p[i-2];
		p[i] = p[i]/(double)i;
	}
}
//===================================================================
//=== Calculate Associated Legendre Polynomial   ====================
//===     only m > 0. In case of -m , multiply (-1)**m*(l-m)!/(l+m)!=
//===                                                               =
//=== Plm = p (l * (l + 1) + m)                                     =
//===================================================================
void alegendre(
	vector<double>& p,
	const double x,
	const int l,
	const int m
)
{
	if(m > l) return;

	p.resize((l + 1) * (l + 1));
	
	double xx = sqrt(1.0-x*x);
	p[0] = 1.0;
	p[3] = xx;
	p[2] = x;
	p[1] = -p[3]*0.5;
	for (int i = 2; i <= l; i++) {
		p[i*(i+1)] = (2.0*(double)i-1.0)*x*p[i*(i-1)]
					-(double)(i-1.0)*p[(i-1)*(i-2)];
		p[i*(i+1)] /= (double)i;
		double tmp = (2.0*(double)i-1.0)*xx;
		int ii = (i == l) ? abs(m) : i;
		for (int j = 1; j <= ii; j++) {
			double tt1 = (j - 1 > i - 1) ? 0.0 : p[i * (i - 1) + j - 1];
			double tt2 = (j > i - 2) ? 0.0 : p[(i - 1) * (i - 2) + j];
			p[i * (i + 1) + j] = tmp * tt1 + tt2;
		}
	}
}
//-----------------------------------------------------------------------
//calculate normalization factor N in the full SH basis formula Y = N * P(cos theta) * sin(m phi) or cos(abs(m) phi)
// N is normalization factor
// P is alegendre
void calcNormalizationCoeff(
	vector<double>& coeff,
	const int degree
)
{
	coeff.clear();
	for(int l=0; l<=degree; l++) {
		for(int m = -l; m<= l; m++) {
			double a = 1;
			for(int k=1; k<=(l-abs(m)); k++) a *= k;
			double b = 1;
			for(int k=1; k<=(l+abs(m)); k++) b *= k;
			coeff.push_back(sqrt(((2.0 * l + 1) * a) / (4.0 * M_PI * b)));
			if(m != 0) coeff.back() *= sqrt(2.0);
		}
	}
}
//-----------------------------------------------------------------------
//calculating SH basis using full formula Y = N * P(cos theta) * sin(m phi) or cos(abs(m) phi)
//N is normalization factor
//P is alegendre
void sphericalHarmonics(
	vector<double>& sh,
	const vector<double>& normalizedCoeffs,
	const int degree,
	const double ct,
	const double phi
)
{
	vector<double> alegendreCoeffs;
    alegendre(alegendreCoeffs, ct, degree, degree); // calculate P | since m = [-l;l] we pass degree to both param
	for (int i = 0; i <= degree; i++) {
		int ii = i*(i+1);
		for(int j=-i; j<=i; j++) {
			int jj = abs(j);
			if(j >= 0) {
                            sh[ii + j] = normalizedCoeffs[ii + j] * alegendreCoeffs[ii + jj] * cos(jj * phi); // coeff is N | tmp is P 
			} else {
                            sh[ii + j] = normalizedCoeffs[ii + j] * alegendreCoeffs[ii + jj] * sin(jj * phi);
			}
		}
	}
}
//-----------------------------------------------------------------------
void sphericalHarmonics(
	vector<double>& sh,
	const int degree,
	const double ct,
	const double phi
)
{
	vector<double> coeff;
	calcNormalizationCoeff(coeff, degree);
	sphericalHarmonics(sh, coeff, degree, ct, phi);
}
//-----------------------------------------------------------------------

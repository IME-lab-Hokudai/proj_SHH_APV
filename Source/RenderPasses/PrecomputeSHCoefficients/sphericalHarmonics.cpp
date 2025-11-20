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
	vector<double>& shBasis,
	const vector<double>& normalizedCoeffs,
	const int degree,
	const double ct,
	const double phi
)
{
	vector<double> alegendreCoeffs;
        alegendre(alegendreCoeffs, ct, degree, degree); // calculate P | since m = [-l;l] we pass degree to both param
	for (int L = 0; L <= degree; L++) {
		int ii = L*(L+1);
		for(int M=-L; M<=L; M++) {
			int absM = abs(M);
			if(M >= 0) {
                            shBasis[ii + M] = normalizedCoeffs[ii + M] * alegendreCoeffs[ii + absM] * cos(absM * phi); // coeff is N | tmp is P 
			} else {
                            shBasis[ii + M] = normalizedCoeffs[ii + M] * alegendreCoeffs[ii + absM] * sin(absM * phi);
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

double CSPhase(int m)
{
    return (std::abs(m) % 2 == 1) ? -1.0 : 1.0;
}

double getNormalizationK(int l, int m)
{
    int m_abs = std::abs(m);

    // Compute factorial terms using std::tgamma (tgamma(n+1) = n!)
    double fact_l_minus_m = std::tgamma(l - m_abs + 1);
    double fact_l_plus_m = std::tgamma(l + m_abs + 1);

    double numerator = (2.0 * l + 1.0) * fact_l_minus_m;
    double denominator = (4.0 * M_PI) * fact_l_plus_m;

    return std::sqrt(numerator / denominator);
}

double get_dP_dx(int l, int m, double x)
{
    int m_abs = std::abs(m);

    // 1. Handle singularities at poles (x -> +/- 1)
    // The term 1/(x^2 - 1) explodes at the poles. We clamp x slightly.
    const double EPSILON = 1e-6;
    if (std::abs(x) > 1.0 - EPSILON)
    {
        x = (x > 0.0) ? (1.0 - EPSILON) : (-1.0 + EPSILON);
    }

    // 2. Evaluate P_l^m(x)
    // C++17 std::assoc_legendre includes the Condon-Shortley phase (-1)^m
    double P_current = std::assoc_legendre(l, m_abs, x);

    // 3. Evaluate P_{l-1}^m(x)
    double P_prev = 0.0;
    if (l - 1 >= m_abs)
    {
        P_prev = std::assoc_legendre(l - 1, m_abs, x);
    }

    // 4. Apply recurrence
    double numerator = ((double)l*x * P_current) - ((double)(l + m_abs) * P_prev);

    double denominator = (x * x) - 1.0;

    // Helper: Condon-Shortley Phase Factor (-1)^m
    // Since std::assoc_legendre OMITS this phase, we must multiply by it
    // if the target basis requires it (which your mismatch suggests).
    double csPhase = CSPhase(m);

    return (numerator / denominator)*csPhase;
}

double computeThetaGradient(int l, int m, double theta, double phi)
{
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    // Calculate normalization K and the Legendre derivative term
    double K = getNormalizationK(l, m);
    double dP_dx_val = get_dP_dx(l, m, cos_theta);

    // Apply specific cases based on the sign of m 
    if (m > 0)
    {
        // Case m > 0: -sqrt(2) * K * cos(m*phi) * sin(theta) * dP/dx
        return -std::sqrt(2.0) * K * std::cos(m * phi) * sin_theta * dP_dx_val;
    }
    else if (m < 0)
    {
        // Case m < 0: -sqrt(2) * K * sin(-m*phi) * sin(theta) * dP/dx
        // Note: paper uses sin(-m * phi)
        return -std::sqrt(2.0) * K * std::sin(-m * phi) * sin_theta * dP_dx_val;
    }
    else
    {
        // Case m = 0: -K * sin(theta) * dP/dx
        return -K * sin_theta * dP_dx_val;
    }
}

double computePhiGradient(int l, int m, double ylmminus)
{
    if (m == 0)
    {
        return 0.0;
    }

    double res = -(double)m * ylmminus;
    return res;
}

//-----------------------------------------------------------------------

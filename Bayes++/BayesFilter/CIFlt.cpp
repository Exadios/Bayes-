/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Intersection Filter.
 */
#include "matSup.h"
#include "CIFlt.h"
#include "models.h"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


CI_filter::CI_filter (size_t x_size, size_t z_initialsize) :
	Extended_filter(x_size),
	S(Empty), SI(Empty)
/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

CI_filter& CI_filter::operator= (const CI_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	return *this;
}


void CI_filter::init ()
{
						// Preconditions
	if (!isPSD (X))
		filter_error ("Xi not PSD");
}

void CI_filter::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 CI_filter::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	RowMatrix temp(f.Fx.size1(), X.size2());
	X = prod_SPD(f.Fx,X, temp) + prod_SPD(f.G, f.q);

	assert_isPSD (X);
	return 1.;
}

void CI_filter::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		S.resize(z_size,z_size);
		SI.resize(z_size,z_size);
	}
}

Bayes_base::Float
 CI_filter::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
/*
 * Iterated Extended Kalman Filter
 * Bar-Shalom and Fortmann p.119 (full scheme)
 * A hard limit is placed on the iterations whatever the
 * the normal terminal condition is to guarantee termination
 * Uncorrelated noise
 */
{
						// ISSUE: Implement simplified uncorrelated noise equations
	size_t z_size = s.size();
	SymMatrix Z(z_size,z_size);
	
	Adapted_Linrz_correlated_observe_model hh(h);
	return observe_innovation (hh, s);
}


Bayes_base::Float
 CI_filter::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* correlated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

	// find omega  //TODO optimise of R.size==1
	SymMatrix invZ(h.Z.size1(),h.Z.size2());
	Float rcond = UdUinversePD (invZ, h.Z);
	rclimit.check_PSD(rcond, "Z not PSD in observe");

	Matrix HTran = trans(h.Hx);
	Matrix HTranInvR = prod(HTran, invZ);
	SymMatrix HTinvRH = prod(HTranInvR, h.Hx);

	SymMatrix invX(X.size1(),X.size2());
	rcond = UdUinversePD (invX, X);
	rclimit.check_PD(rcond, "X not PD in observe");

std::cout << X <<std::endl;
std::cout << invX <<std::endl;
std::cout << HTinvRH <<std::endl;

	double omega = Omega(invX, HTinvRH, X);
std::cout << omega <<std::endl;

	/* calculate predicted innovation */
	Matrix PHTran = prod(X, HTran);
	Matrix HPHTran = prod(h.Hx, PHTran);
std::cout << HPHTran <<std::endl;
	Matrix oHPHTran = HPHTran * (1.-omega);
	SymMatrix oR = h.Z * omega;

	S = oHPHTran + oR;
std::cout << S <<std::endl;

	/* test for fault ??*/
	//    error = faultTest( S, obsErr, obsNum  );

	// Inverse innovation covariance
	rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

	Matrix oPHTran = PHTran*(1.0-omega);
	Matrix K = prod(oPHTran, SI);
std::cout << K <<std::endl;

	// State update
	x += prod(K, s);
	// Invserse covariance
	SymMatrix oInvCov = invX * omega;
	SymMatrix oHTinvRH = HTinvRH*(1.0-omega);
	invX = oInvCov +oHTinvRH;
	// Covariance
	rcond = UdUinversePD (X, invX);
	rclimit.check_PD(rcond, "inverse covariance not PD in observe");
	return rcond;
}


CI_filter::Float CI_filter::Omega(const SymMatrix& Ai_ext, const SymMatrix& Bi_ext, const SymMatrix& A_ext)
/*
 * Compute CI Omega value
 *  Iterative optimization algorithm from the authors of reference [1]
 */
{
	using std::sqrt;
	const int MAX_STATES=10;	// Nast fixed value for C implementation
	const double MINDOUBLE=1E-38;
	int n = A_ext.size1();
	int i, j, k, m, count;
	double a=0.0, b, c, f, g, p, q, r, s;
	double d[MAX_STATES], e[MAX_STATES], t=0.0, dt=0.0, X[MAX_STATES][MAX_STATES];

	Matrix Ai = Ai_ext;			// Copy extern matrices internal work matrices
	Matrix Bi = Bi_ext;
	Matrix A  = A_ext;

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			for (k=0,s=0; k<j; k++) s += X[i][k] * X[j][k];
			if (j<i) X[i][j] = (A[i][j] - s)/X[j][j];
			else X[i][j] = sqrt(fabs(A[i][i] - s));
		}
	}
	for (i=0; i<n; i++) for (j=i; j<n; j++) A[i][j] = X[j][i];
	for (i=0; i<n; i++) {
		for (j=i; j<n; j++) e[j] = A[i][j];
		for (j=0; j<n; j++) {
			for (k=i,s=0.0; k<n; k++) s += e[k]*Bi[k][j];
			A[i][j] = s;
		}
	}
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) e[j] = A[i][j];
		for (j=0; j<n; j++) {
			for (k=j,s=0.0; k<n; k++) s += e[k]*X[k][j];
			A[i][j] = s;
		}
	}
	for (i=(n-1); i>0; i--) {
		m = i - 1;  r = s = 0.0;
		if (m > 0) {
			for (k=0; k<=m; k++) s += fabs(A[i][k]);
			if (s == 0.0) e[i] = A[i][m];
			else {
				for (k=0; k<=m; k++) {
					A[i][k] /= s;
					r += A[i][k]*A[i][k];
				}
				g = ((f = A[i][m]) >= 0.0 ? -sqrt(r) : sqrt(r));
				e[i] = s * g;  r -= f*g;
				A[i][m] = f - g;  f = 0.0;
				for (j=0; j<=m; j++) {
					for (k=0,g=0.0; k<=j; k++) g += A[j][k]*A[i][k];
					for (k=j+1; k<=m; k++) g += A[k][j]*A[i][k];
					e[j] = g/r;
					f += e[j]*A[i][j];
				}
				q = f/(r+r);
				for (j=0; j<=m; j++) {
					f = A[i][j];
					g = (e[j] -= q*f);
					for (k=0; k<=j; k++) A[j][k]-=(f*e[k]+g*A[i][k]);
				}
			}
		}
		else e[i] = A[i][m];
		d[i] = r;
	}
	for (i=0; i<n; i++) d[i] = A[i][i];
	for (i=1; i<n; i++) e[i-1] = e[i];
	e[n-1] = 0.0;
	for (j=0; j<n; j++) {
		count = 0;
		do {
			for (m=j; m<n-1; m++) {
				q = fabs(d[m]) + fabs(d[m+1]);
				if ((double)(fabs(e[m])+q) == q) break;
			}
			if (m != j) {
				if (count++ > sizeof(double)*4) goto Iloop;
				g = (d[j+1]-d[j])/(2.0*e[j]);
				r = sqrt(g*g+1.0);
				g = d[m]-d[j]+e[j]/(g+((g<0)?-fabs(r):fabs(r)));
				s = c = 1.0; p = 0.0;
				for (i=m-1; i>=j; i--) {
					f = s*e[i]; b = c*e[i];
					e[i+1] = (r = sqrt(f*f+g*g));
					if (r == 0.0) {
						d[i+1] -= p; e[m] = 0.0;
						break;
					}
					s = f/r; c = g/r;
					g = d[i+1] - p;
					r = (d[i]-g)*s + 2.0*c*b;
					d[i+1] = g + (p=s*r);
					g = c*r - b;
				}
				if (r == 0.0 && i >= 0) continue;
				d[j] -= p; e[j] = g; e[m] = 0.0;
			}
		} while (m != j);
	}
Iloop:
	for (i=0,b=c=0.0; i<n; i++) {
		if (d[i]==0.0) j = 1;
		else b += (1-d[i])/d[i];
		c += (1-d[i]);
	}
	if (j) b = fabs(c)+2*(a=sqrt(MINDOUBLE));
	if (fabs(b-c)<a)
		return(1.0);
	if (b*c>0.0) {
		if (b<0.0)
			return(0.0);
		else
			return(1.0);
	}
	if (b > 0.0) {c = 0.0; b = 1.0;}
	else {b = 0.0; c = 1.0;}
	for (i=0; i<n; i++) {
		t += (q = 2*(1-d[i])/(1+d[i]));
		dt -= q*q;
	}
	p = q = 1.0; g = 0.5;
	for (j=0; j<(32*sizeof(double)); j++) {
		if ((((g-c)*dt-t)*((g-b)*dt-t) >= 0.0)
			|| (fabs(2.0*t) > fabs(q*dt))) {
				q = p;
				p = 0.5*(c-b);
				g = b + p;
				if (b == g)
					goto Iwrap;
			}
		else {
			q = p; p = t/dt; q = g; g -= p;
			if (q == g)
				goto Iwrap;
		}
		if (fabs(p) < a)
			goto Iwrap;
		t = dt = 0.0;
		for (i=0; i<n; i++) {
			t += (q = (1-d[i])/(g*(1-d[i])+d[i]));
			dt -= q*q;
		}
		if (t < 0.0) b = g;
		else c = g;
	}
Iwrap:
	if (1.0-g<=2.*a) return(1.0);
	if (g<=2.*a) return(0.0);
	return(g);
}

}//namespace

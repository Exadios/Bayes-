/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Root Filter.
 */
#include "bayesFlt.hpp"
#include "infRtFlt.hpp"
#include "matSup.hpp"
#include "uLAPACK.hpp"	// Common LAPACK interface


/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;

							// Simple type to get upper triangular part of of a matrix
typedef ublas::triangular_adaptor<const Matrix, ublas::upper> UpperTri;


Information_root_filter::Information_root_filter (size_t x_size, size_t /*z_initialsize*/) :
		Extended_filter(x_size),
		r(x_size), R(x_size,x_size)
/*
 * Set the size of things we know about
 */
{
	// No temporaries to make use of z_initialsize
}

Information_root_info_filter::Information_root_info_filter (size_t x_size, size_t z_initialsize) :
		Information_root_filter (x_size, z_initialsize),
		y(x_size), Y(x_size,x_size)
{}

void Information_root_filter::init ()
/*
 * Inialise the filter from x,X
 * Predcondition:
 *		X is PSD
 * Postcondition:
 *		inv(R)*inv(R)' = X
 *		r = R*x
 */
{
						// Information Root
	Float rcond = UCfactor (R, X);
	rclimit.check_PSD(rcond, "Xi not PSD");
	(void)UTinverse(R);
						// Information Root state r=R*x
	r.assign (prod(R,x));
}

void Information_root_filter::init_information (const Vec& yi, const SymMatrix& Yi)
/*
 * Special Initialisation directly from Information form
 * Postcondition:
 *		R'*R = Yi
 *		r = inv(R)'*y
 */
{
					// Temporary R Matrix for factorisation
	const size_t n = x.size();
	LTriMatrix LC(n,n);
					// Information Root
	Float rcond = LdLfactor (LC, Yi);
	rclimit.check_PSD(rcond, "Yi not PSD");

	{				// Lower triangular Choleksy factor of LdL'
		size_t i,j;
		for (i = 0; i < n; ++i)
		{
			using namespace std;		// for sqrt
			LTriMatrix::value_type sd = LC(i,i);
			sd = sqrt(sd);
			LC(i,i) = sd;
						// Multiply columns by square of non zero diagonal
			for (j = i+1; j < n; ++j)
			{
				LC(j,i) *= sd;		// TODO Use Row operation
			}
		}
	}
	R = FM::trans(LC);			// R = (L*sqrt(d))'

	UTriMatrix RI(n,n);
	RI = R;
	(void)UTinverse(RI);
	r.assign (prod(FM::trans(RI),yi));
}

void Information_root_filter::update ()
/*
 * Recompute x,X from r,R
 * Precondition:
 *		r(k|k),R(k|k)
 * Postcondition:
 *		r(k|k),R(k|k)
 *		x = inv(R)*r
 *		X = inv(R)*inv(R)'
 */
{
	UTriMatrix RI(x.size(), x.size());
	RI = R;			// Invert Cholesky factor
	bool singular = UTinverse(RI);
	if (singular)
		filter_error ("R not PD");

	X.clear();
	mult_SPDi(RI, X);		// X = RI*RI'
	x.assign (prod(RI,r));

	assert_isPSD (X);
}

void Information_root_info_filter::update ()
/*
 * Recompute x,X y,Y from r,R
 * Precondition:
 *		r(k|k),R(k|k)
 * Postcondition:
 *		r(k|k),R(k|k)
 *		x = inv(R)*r
 *		X = inv(R)*inv(R)'
 *		Y = R'*R
 *		y = Y*x
 */
{
	Information_root_filter::update();
	Y.clear();
	mult_SPDTi(R, Y);
	y.assign (prod(Y,x));
}


Bayes_base::Float
 Information_root_filter::predict (Linrz_predict_model& f)
/* Linrz Prediction: goes back to state to compute
 * Precondition:
 *		r(k|k),R(k|k)
 * Postcondition:
 *		r(k+1|k) = R(k+1|k) * f(x(k|k)
 *		R(k+1|k) computed using QR decomposition after inverse(Fx) is computed
 *
 * Uses LAPACK geqrf for QR decomposition (without PIVOTING)
 */
{
	const size_t q_size = f.q.size();
	const size_t x_size = x.size();

	update ();			// x is required for f(x);

						// Require inverse(Fx)
	ColMatrix FxI(x_size, x_size);
	{								// LU factorise with pivots
		ColMatrix FxLU(x_size, x_size);
		FxLU = f.Fx;
		LAPACK::pivot_t ipivot(x_size);		// Pivots initialised to zero
		for (LAPACK::pivot_t::iterator i = ipivot.begin(); i != ipivot.end(); ++i)
			*i = 0;
		int info = LAPACK::getrf(FxLU, ipivot);
		if (info < 0)
			filter_error ("Fx not LU factorisable");

		FM::identity(FxI);				// Invert
		info = LAPACK::getrs('N', FxLU, ipivot, FxI); 
		if (info != 0)
			filter_error ("Predict Fx not LU invertable");
	}
						// Require Root of correlated predict noise (may be semidefinite)
	Matrix Qr(q_size, q_size);		// ISSUE Requires a diagonal matrix
	Qr.clear();
	for (Vec::const_iterator i = f.q.begin(); i != f.q.end(); ++i)
	{
		using namespace std;
		size_t ii = i.index();
		if (*i < 0.)
			filter_error ("Negative q");
		Qr(ii,ii) = sqrt(*i);
	}
						// Form Augmented matrix for factorisation
	ColMatrix A(x_size*2, x_size*2);		// Column major required for LAPACK, also this property is using in indexing
	FM::identity (A);						// Prefill with identity for topleft and zero's in off diagonals

	A.sub_matrix(x_size,x_size*2, 0,x_size).assign (prod(R, prod(FxI, prod(f.G,Qr))));
	A.sub_matrix(x_size,x_size*2, x_size,x_size*2).assign (prod(R, FxI));

						// Calculate factorisation so we have and upper triangular R
	Vec tau(x_size*2);
	int info = LAPACK::geqrf (A, tau);
	if (info != 0)
			filter_error ("Predict no QR factor");
						// Extract the roots, junk in strict lower triangle
	R = UpperTri(A.sub_matrix(x_size,x_size*2, x_size,x_size*2));
					
	r.assign (prod(R,f.f(x)));	// compute r using f(x)

	return UCrcond(R);	// compute rcond of result
}


Bayes_base::Float
 Information_root_filter::predict (Linear_invertable_predict_model& f)
/* Linear Prediction using linear invertable model
 * Precondition:
 *		r(k|k),R(k|k)
 * Postcondition:
 *		r(k+1|k) = R(k+1|k) * f(x(k|k)
 *		R(k+1|k) computed using QR decomposition using inverse(Fx)
 *
 * Uses LAPACK geqrf for QR decomposition (without PIVOTING)
 */
{
	const size_t q_size = f.q.size();
	const size_t x_size = x.size();

						// Require Root of correlated predict noise (may be semidefinite)
	DiagMatrix Qr(q_size, q_size);
	Qr.clear();
	for (Vec::const_iterator i = f.q.begin(); i != f.q.end(); ++i)
	{
		using namespace std;
		size_t ii = i.index();
		if (*i < 0.)
			filter_error ("Negative q");
		Qr(ii,ii) = sqrt(*i);
	}
						// Form Augmented matrix for factorisation
	ColMatrix A(x_size*2, x_size*2+1);		// Column major required for LAPACK, also this property is using in indexing
	FM::identity (A);						// Prefill with identity for topleft and zero's in off diagonals

	A.sub_matrix(x_size,x_size*2, 0,x_size).assign (prod(R, prod(f.inv.Fx, prod(f.G,Qr))));
	A.sub_matrix(x_size,x_size*2, x_size,x_size*2).assign (prod(R, f.inv.Fx));
	A.sub_column(x_size,x_size*2, x_size*2).assign (r);

						// Calculate factorisation so we have and upper triangular R
	Vec tau(x_size*2);
	int info = LAPACK::geqrf (A, tau);
	if (info != 0)
			filter_error ("Predict no QR factor");

						// Extract the roots, junk in strict lower triangle
	R = UpperTri(A.sub_matrix(x_size,x_size*2, x_size,x_size*2));
	r = A.sub_column(x_size,x_size*2, x_size*2);

	return UCrcond(R);	// compute rcond of result
}


Bayes_base::Float Information_root_filter::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* Extended linrz correlated observe
 * Precondition:
 *		r(k+1|k),R(k+1|k)
 * Postcondition:
 *		r(k+1|k+1),R(k+1|k+1)
 *
 * Uses LAPACK geqrf for QR decomposigtion (without PIVOTING)
 * ISSUE correctness of linrz form needs validation
 */
{
	const size_t x_size = x.size();
	const size_t z_size = s.size();
						// Size consistency, z to model
	if (z_size != h.Z.size1())
		filter_error("observation and model size inconsistent");

						// Require Inverse of Root of uncorrelated observe noise
	UTriMatrix Zir(z_size,z_size);
	Float rcond = UCfactor (Zir, h.Z);
	rclimit.check_PSD(rcond, "Z not PSD");
	UTinverse(Zir);
						// Form Augmented matrix for factorisation
	ColMatrix A(x_size+z_size, x_size+1);	// Column major required for LAPACK, also this property is using in indexing
	A.sub_matrix(0,x_size, 0,x_size).assign (R);
	A.sub_matrix(x_size,x_size+z_size, 0,x_size).assign (prod(Zir, h.Hx));
	A.sub_column(0,x_size, x_size).assign (r);
	A.sub_column(x_size,x_size+z_size, x_size).assign (prod(Zir, s+prod(h.Hx,x)));

						// Calculate factorisation so we have and upper triangular R
	Vec tau(x_size+1);
	int info = LAPACK::geqrf (A, tau);
	if (info != 0)
			filter_error ("Observe no QR factor");
						// Extract the roots, junk in strict lower triangle
	R = UpperTri(A.sub_matrix(0,x_size, 0,x_size));
	r = A.sub_column(0,x_size, x_size);

	return UCrcond(R);	// compute rcond of result
}


Bayes_base::Float Information_root_filter::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
/* Extended linrz uncorrelated observe
 * Precondition:
 *		r(k+1|k),R(k+1|k)
 * Postcondition:
 *		r(k+1|k+1),R(k+1|k+1)
 *
 * Uses LAPACK geqrf for QR decomposigtion (without PIVOTING)
 * ISSUE correctness of linrz form needs validation
 * ISSUE Efficiency. Product of Zir can be simplified
 */
{
	const size_t x_size = x.size();
	const size_t z_size = s.size();
						// Size consistency, z to model
	if (z_size != h.Zv.size())
		filter_error("observation and model size inconsistent");

						// Require Inverse of Root of uncorrelated observe noise
	DiagMatrix Zir(z_size,z_size);
	Zir.clear();
	for (size_t i = 0; i < z_size; ++i)
	{
		using namespace std;
		Zir(i,i) = 1./ sqrt(h.Zv[i]);
	}
						// Form Augmented matrix for factorisation
	ColMatrix A(x_size+z_size, x_size+1);	// Column major required for LAPACK, also this property is using in indexing
	A.sub_matrix(0,x_size, 0,x_size).assign (R);
	A.sub_matrix(x_size,x_size+z_size, 0,x_size).assign (prod(Zir, h.Hx));
	A.sub_column(0,x_size, x_size).assign (r);
	A.sub_column(x_size,x_size+z_size, x_size).assign (prod(Zir, s+prod(h.Hx,x)));

						// Calculate factorisation so we have and upper triangular R
	Vec tau(x_size+1);
	int info = LAPACK::geqrf (A, tau);
	if (info != 0)
			filter_error ("Observe no QR factor");
						// Extract the roots, junk in strict lower triangle
	R = UpperTri(A.sub_matrix(0,x_size, 0,x_size));
	r = A.sub_column(0,x_size, x_size);

	return UCrcond(R);	// compute rcond of result
}

}//namespace

#ifndef _BAYES_FILTER_MATRIX_SUPPORT
#define _BAYES_FILTER_MATRIX_SUPPORT

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Matrix support functions for filter classes
 *  Members of the Bayesian_filter_matrix namespace are used in the
 *  interface of Bayesian_filter and for internal operations.
 *  Be aware! These functions and their implemenation are more likely to change
 *   then those in Bayesian_filter.
 */

#include "matSupSub.h"		// Expect to find the actual matrix support headers elsewhere

/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{


/*
 * Assertion support
 */
#ifndef NDEBUG
void assert_isSymmetric (const SymMatrix &M);
void assert_isPSD (const SymMatrix &M);
#else
inline void assert_isSymmetric (const SymMatrix &M) {}
inline void assert_isPSD (const SymMatrix &M) {}
#endif

/*
 * Local support functions
 */
bool isPSD (const SymMatrix &M);
bool isSymmetric (const SymMatrix &M);
void forceSymmetric (SymMatrix &M, bool bUpperToLower = false);

/*
 * UdU' and LdL' and UU' Cholesky Factorisation and function
 * Very important to manipulate PD and PSD matrices
 */
RowMatrix::value_type UdUrcond (const RowMatrix& UD);
SymMatrix::value_type UdUrcond (const SymMatrix& UD);
UTriMatrix::value_type UCrcond (const UTriMatrix& U);
Vec::value_type UdUrcond_vec (const Vec& D);
SymMatrix::value_type UdUdet (const SymMatrix& UD);

RowMatrix::value_type UdUfactor_varient1 (RowMatrix& M, size_t n);
RowMatrix::value_type UdUfactor_varient2 (RowMatrix& M, size_t n);
inline RowMatrix::value_type UdUfactor (RowMatrix& M, size_t n)
{	return UdUfactor_varient2(M,n);
}
LTriMatrix::value_type LdLfactor (LTriMatrix& M, size_t n);
UTriMatrix::value_type UCfactor (UTriMatrix& M, size_t n);

RowMatrix::value_type UdUfactor (RowMatrix& UD, const SymMatrix& M);
LTriMatrix::value_type LdLfactor (LTriMatrix& LD, const SymMatrix& M);
UTriMatrix::value_type UCfactor (UTriMatrix& U, const SymMatrix& M);

bool UdUinverse (RowMatrix& UD);
bool UTinverse (UTriMatrix& U);
void UdUrecompose_transpose (RowMatrix& M);
void UdUrecompose (RowMatrix& M);
void UdUfromUCholesky (RowMatrix& U);
void UdUseperate (RowMatrix& U, Vec& d, const RowMatrix& UD);
void Lzero (RowMatrix& M);
void Uzero (RowMatrix& M);

/*
 * Functions using UdU factorisation:
 *  inverse of Positive Definate matrix
 */
SymMatrix::value_type UdUinversePD (SymMatrix& M);
SymMatrix::value_type UdUinversePD (SymMatrix& M, SymMatrix::value_type& detM);
SymMatrix::value_type UdUinversePD (SymMatrix& MI, const SymMatrix& M);
SymMatrix::value_type UdUinversePD (SymMatrix& MI, SymMatrix::value_type& detM, const SymMatrix& M);


/*
 * Matrix Support Operations
 */
template <class BaseL, class BaseR>
void subcopy (const FMMatrix<BaseL>& L, FMMatrix<BaseR>& R)
{	// Copy L into top left block of R
	using namespace ublas;
	matrix_range<BaseR> Rtopleft(R, range(0,L.size1()), range(0,L.size2()));
	Rtopleft = L;
}

template <class Base>
void identity(FMMatrix<Base>& I)
// Fill as identity matrix
{
	I.clear();
							// Set common diagonal elements
	size_t common_size = std::min(I.size1(),I.size2());
	typedef typename Base::value_type Base_value_type;
	ublas::matrix_vector_range<Base>(I, ublas::range(0,common_size), ublas::range(0,common_size)) = 
		ublas::scalar_vector<Base_value_type>(common_size, 1.);
}


/*
 * Symmetric Positive (Semi) Definate multiply:	X*S*X'
 *  Couldn't find a better name for this operation. Result is actually only PD if S is
 */

inline Vec::value_type mult_SPD (const Vec& x, const SymMatrix& S)
/*
 * Symmetric Positive (Semi) Definate multiply:	p = x*S*x'
 * Optimised to exploit the symmetry of the computation of x * x' and result p
 * S must be Symetric (SPD -> p is SPD)
 */
{
	Vec::value_type p = 0.;
	Vec::const_iterator xi = x.begin();
	for (SymMatrix::const_iterator1 Si = S.begin1(); Si != S.end1(); ++Si, ++xi)
	{
		p += *xi * ublas::inner_prod(SymMatrix::row(Si), x);
	}
	return p;
}

inline Vec::value_type mult_SPD (const Vec& x, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate multiply:	p = x*diag_matrix(s)*x'
 * Optimised to exploit the symmetry of the computation of x * x' and result p
 */
{
	Vec::value_type p = 0.;
	Vec::const_iterator xi = x.begin();
	for (Vec::const_iterator si = s.begin(); si != s.end(); ++si, ++xi)
	{
		p += (*xi)*(*xi) * (*si);
	}
	return p;
}

template <class MatrixX>
void mult_SPD (const MatrixX& X, const Vec& s, SymMatrix& P)
/*
 * Symmetric Positive (Semi) Definate multiply:	P += X*diag_matrix(s)*X'
 * Optimised to exploit the symmetry of the computation of X * X' and symmetry of result P
 */
{
	Vec::const_iterator si, send = s.end();
	SymMatrix::iterator1 Pa = P.begin1();
	typename MatrixX::const_iterator1 Xa = X.begin1();
	const typename MatrixX::const_iterator1 Xend = X.end1();
	typename MatrixX::const_iterator1 Xb;

	// P(a,b) = X.row(a) * X.row(b)
	for (; Xa != Xend; ++Xa,++Pa)			// Iterate Rows
	{		
		Xb = Xa;							// Start at the row Xa only one triangle of symetric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb != Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab;				// Simple vector operation
			for (si = s.begin(); si != send; ++si) {
				Vec::size_type i = si.index();
				p += *(Xa.begin()+i) * (*si) * *(Xb.begin()+i);	// TODO: use iterator on *Xa and *Xb
			}
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}

template <class MatrixX>
void mult_SPDi (const MatrixX& X, SymMatrix& P)
/*
 * Symmetric Positive (Semi) Definate multiply:	P += X*X'
 *  Result is always PD
 * Optimised to exploit the symmetry of the computation of X * X' and symmetry of result P
 */
{
	SymMatrix::iterator1 Pa = P.begin1();
	typename MatrixX::const_iterator1 Xa = X.begin1();
	const typename MatrixX::const_iterator1 Xend = X.end1();
	typename MatrixX::const_iterator1 Xb;

	// P(a,b) = X.row(a) * X.row(b)
	for (; Xa != Xend; ++Xa,++Pa)			// Iterate vectors
	{		
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb != Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab;				// Simply multiple row Xa by row Xb
			p += ublas::inner_prod( MatrixX::row(Xa), MatrixX::row(Xb) );
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}

template <class MatrixX>
void mult_SPD (const MatrixX& X, const SymMatrix& S, SymMatrix& P, Vec& stemp)
/*
 * Symmetric Positive (Semi) Definate multiply: P += X*S*X'
 *  A temporary Vec is required of same size as S
 * Optimised to exploit the symmetry of the computation of X * X' and symmetry of result P
 * S must be Symetric (SPD -> P remains SPD)
 */
{
	SymMatrix::iterator1 Pa = P.begin1();
	typename MatrixX::const_iterator1 Xa = X.begin1();
	const typename MatrixX::const_iterator1 Xend = X.end1();
	typename MatrixX::const_iterator1 Xb;

	// P(a,b) = X.row(a) * S * X'.col(b) = X.row(a) * S * X.row(b)
	//        = (S * X.row(a)')' * X.row(b)

	for (; Xa != Xend; ++Xa,++Pa)			// Iterate Rows
	{
		stemp.assign( ublas::prod(S, MatrixX::row(Xa)) );			// treated as a column vector
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb!= Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab + ublas::inner_prod(stemp, MatrixX::row(Xb));
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}


template <class MatrixX>
void mult_SPDT (const MatrixX& X, const SymMatrix& S, SymMatrix& P, Vec& stemp)
/*
 * Symmetric Positive (Semi) Definate multiply: P += X'*S*X
 *  A temporary Vec is required of same size as S
 * Optimised to exploit the symmetry of the computation of X' * X and symmetry of result P
 * S must be Symetric (SPD -> P remains SPD)
 */
{
	SymMatrix::iterator1 Pa = P.begin1();
	typename MatrixX::const_iterator2 Xa = X.begin2();
	const typename MatrixX::const_iterator2 Xend = X.end2();
	typename MatrixX::const_iterator2 Xb;

	// P(a,b) = X'.row(a) * S * X.col(b) = X.col(a) * S * X.col(b)
	//        = X.col(b) * S * X.col(a)

	for (; Xa != Xend; ++Xa,++Pa)			// Iterate Rows
	{
		stemp = ublas::prod(S, MatrixX::column(Xa));			// treated as a column vector
		Xb = Xa;							// Start at the column Xa, only one triangle of symertric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb != Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab + ublas::inner_prod(stemp, MatrixX::column(Xb));
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}

template <class MatrixX>
void mult_SPDT (const MatrixX& X, const Vec& s, SymMatrix& P)
/*
 * Symmetric Positive (Semi) Definate multiply:	P += X'*diag_matrix(s)*X
 * Optimised to exploit the symmetry of the computation of X' * X and symmetry of result P
 */
{
	Vec::const_iterator si, send = s.end();
	SymMatrix::iterator1 Pa = P.begin1();
	typename MatrixX::const_iterator2 Xa = X.begin2();
	const typename MatrixX::const_iterator2 Xend = X.end2();
	typename MatrixX::const_iterator2 Xb;

	// P(a,b) = X.col(a) * X.col(b)
	for (; Xa != Xend; ++Xa,++Pa)			// Iterate vectors
	{		
		Xb = Xa;							// Start at the row Xa only one triangle of symertric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb != Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab;				// Simple vector operation
			for (si = s.begin(); si != send; ++si) {
				Vec::size_type i = si.index();
				p += *(Xa.begin()+i) * (*si) * *(Xb.begin()+i);	// TODO: use iterator on *Xa and *Xb
			}
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}

template <class MatrixX>
void mult_SPDTi (const MatrixX& X, SymMatrix& P)
/*
 * Symmetric Positive (Semi) Definate multiply:	P += X'*X
 *  Result is always PD
 * Optimised to exploit the symmetry of the computation of X' * X and symmetry of result P
 */
{
	SymMatrix::iterator1 Pa = P.begin1();
	typename MatrixX::const_iterator2 Xa = X.begin2();
	const typename MatrixX::const_iterator2 Xend = X.end2();
	typename MatrixX::const_iterator2 Xb;

	// P(a,b) = X.row(a) * X.row(b)
	for (; Xa != Xend; ++Xa,++Pa)			// Iterate vectors
	{		
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb != Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab;				// Simply multiple col Xa by col Xb
			p += ublas::inner_prod( MatrixX::column(Xa), MatrixX::column(Xb) );
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}


/*
 * Operations that create temporaies
 */
#if !defined(BAYESFILTER_NOTEMPOPS)


inline SymMatrix mult_SPD (const RowMatrix& X, const SymMatrix& S)
/*
 * Symetric Positive (Semi) Definate multiply: return X*S*X'
 */
{
	SymMatrix P(X.size1(),X.size1());
	P.clear();
	Vec stemp(S.size1());
	mult_SPD(X,S, P, stemp);
	return P;
}

inline SymMatrix mult_SPD (const RowMatrix& X, const Vec& s)
/*
 * Symetric Positive (Semi) Definate multiply: return X*diag_matrix(s)*X'
 */
{
	SymMatrix P(X.size1(),X.size1());
	P.clear();
	mult_SPD(X,s, P);
	return P;
}

inline SymMatrix mult_SPD (const RowMatrix& X)
/*
 * Symetric Positive (Semi) Definate multiply:	return X*X'
 */
{
	SymMatrix P(X.size1(),X.size1());
	P.clear();
	mult_SPD(X, P);
	return P;
}

inline SymMatrix mult_SPDT (const RowMatrix& X, const SymMatrix& S)
/*
 * Symetric Positive (Semi) Definate multiply: return X'*S*X
 */
{
	SymMatrix P(X.size2(),X.size2());
	P.clear();
	Vec stemp(S.size1());
	mult_SPDT(X,S, P, stemp);
	return P;
}

inline SymMatrix mult_SPDT (const RowMatrix& X, const Vec& s)
/*
 * Symetric Positive (Semi) Definate multiply:	return X'*diag_matrix(s)*X
 */
{
	SymMatrix P(X.size2(),X.size2());
	P.clear();
	mult_SPDT(X,s, P);
	return P;
}

inline SymMatrix mult_SPDT (const RowMatrix& X)
/*
 * Symetric Positive (Semi) Definate multiply:	return X'*X
 */
{
	SymMatrix P(X.size2(),X.size2());
	P.clear();
	mult_SPDT(X, P);
	return P;
}

#endif

}//namespace

#endif

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
 *
 * Return values:
 *  Many algorithms return a value_type which is a reciprocal condition number
 *  These values are documented for each algorithm and are important way to
 *  determine the validity of the results
 */
RowMatrix::value_type UdUrcond (const RowMatrix& UD);
SymMatrix::value_type UdUrcond (const SymMatrix& UD);
UTriMatrix::value_type UCrcond (const UTriMatrix& U);
Vec::value_type UdUrcond_vec (const Vec& D);
SymMatrix::value_type UdUdet (const SymMatrix& UD);

// In-place factorisations
RowMatrix::value_type UdUfactor_variant1 (RowMatrix& M, size_t n);
RowMatrix::value_type UdUfactor_variant2 (RowMatrix& M, size_t n);
inline RowMatrix::value_type UdUfactor (RowMatrix& M, size_t n)
{	return UdUfactor_variant2(M,n);
}
LTriMatrix::value_type LdLfactor (LTriMatrix& M, size_t n);
UTriMatrix::value_type UCfactor (UTriMatrix& M, size_t n);

// Copy factorisations
RowMatrix::value_type UdUfactor (RowMatrix& UD, const SymMatrix& M);
LTriMatrix::value_type LdLfactor (LTriMatrix& LD, const SymMatrix& M);
UTriMatrix::value_type UCfactor (UTriMatrix& U, const SymMatrix& M);

// Factor manipulations
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
 *  inverse of Positive Definate matrix retruning rcond
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
 *  Optimised product for dense matrices
*   Exploits the symmetry of the computation of X * X' and symmetry of result
 *  Couldn't find a better name for this operation. Result is actually only PD if S is PD
 */

template <class X, class S>
struct mult_SPD_matrix_matrix_traits
{	// Provide ET result type for prod(X,prod(S,trans(X))
	typedef ublas::matrix_unary2_traits<X, ublas::scalar_identity<typename X::value_type> >::result_type  XT_type;
	typedef ublas::matrix_matrix_binary_traits<typename S::value_type, S, 
                                        XT_type::value_type, XT_type>::result_type  SXT_type;
	typedef ublas::matrix_matrix_binary_traits<typename X::value_type, X, 
                                        SXT_type::value_type, SXT_type>::result_type  XSXT_type;
// ISSUE Final part of mult_SPD ET is symmetric. How can this be done?
// For now skip this
	typedef XSXT_type  result_type;
};

/* VERY BAD. The temporay created by prod<T> is reference and then goes out of scope
template <class X, class T, class S>
struct mult_SPD_matrix_temp_matrix_traits
{	// Provide ET result type for prod(X,prod<T>(S,trans(X))
	typedef ublas::matrix_unary2_traits<X, ublas::scalar_identity<typename X::value_type> >::result_type  XT_type;
	typedef T  SXT_type;
	typedef ublas::matrix_matrix_binary_traits<typename X::value_type, X, 
                                        typename SXT_type::value_type, SXT_type>::result_type  XSXT_type;
// ISSUE Final part of mult_SPD ET is symmetric. How can this be done?
// For now skip this
	typedef XSXT_type  result_type;
};
*/ 
/*template <>
struct mult_SPD_traits<RowMatrix,Vec>
{	// Provide return value ET type for mult_SPD
	typedef RowMatrix X; typedef Vec S
};*/

inline Vec::value_type mult_SPD (const Vec& x, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate multiply:	p = x*diag_matrix(s)*x'
 * Optimised to exploit the symmetry of the computation of x * x' and result p
 * ISSUE: not sparse efficient
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
			p += inner_prod( MatrixX::rowi(Xa), MatrixX::rowi(Xb) );
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}

template <class MatrixX>
void mult_SPD (const MatrixX& X, const SymMatrix& S, SymMatrix& P, Vec& stemp)
/*
 * Symmetric Positive (Semi) Definate multiply: P += X*S*X'
 *  A temporary Vec is required of same size as S
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
		stemp.assign( ublas::prod(S, MatrixX::rowi(Xa)) );			// treated as a column vector
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb!= Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab + inner_prod(stemp, MatrixX::rowi(Xb));
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}


inline Vec::value_type mult_SPDT (const Vec& x, const SymMatrix& S)
/*
 * Symmetric Positive (Semi) Definate multiply:	p = x'*S*x
 */
{
	return inner_prod (x, ublas::prod(S,x) );
}

template <class MatrixX>
void mult_SPDT (const MatrixX& X, const SymMatrix& S, SymMatrix& P, Vec& stemp)
/*
 * Symmetric Positive (Semi) Definate multiply: P += X'*S*X
 *  A temporary Vec is required of same size as S
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
		stemp = ublas::prod(S, MatrixX::columni(Xa));			// treated as a column vector
		Xb = Xa;							// Start at the column Xa, only one triangle of symertric result required
		SymMatrix::iterator2 Pab= Pa.begin();
		for (; Xb != Xend; ++Xb,++Pab)
		{
			SymMatrix::value_type p = *Pab + inner_prod(stemp, MatrixX::columni(Xb));
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}

template <class MatrixX>
void mult_SPDT (const MatrixX& X, const Vec& s, SymMatrix& P)
/*
 * Symmetric Positive (Semi) Definate multiply:	P += X'*diag_matrix(s)*X
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
			p += inner_prod( MatrixX::columni(Xa), MatrixX::columni(Xb) );
			*Pab = p; P(Pab.index2(),Pab.index1()) = p;
		}
	}
}


/*
 * Operations that create temporaies
 */
#if !defined(BAYESFILTER_NOTEMPOPS)




inline
mult_SPD_matrix_matrix_traits<RowMatrix,SymMatrix>::result_type 
	mult_SPD (const RowMatrix& X, const SymMatrix& S)
/*
 * Symmetric Positive (Semi) Definate multiply: return X*S*X'
 */
{
	return prod( X, prod(S,trans(X)) );
// ISSUE: need trait for 
//	symmetric_adaptor<XSXT_type::matrix_type> prod (X, prod(S,trans(X)));
/*	SymMatrix P(X.size1(),X.size1());
	P.clear();
	Vec stemp(S.size1());
	mult_SPD(X,S, P, stemp);
	return P;*/
}

inline SymMatrix mult_SPD (const RowMatrix& X, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate multiply: return X*diag_matrix(s)*X'
 */
{
//    return prod (X, prod(s,trans(X))); Need a solution to create diag matrix
	SymMatrix P(X.size1(),X.size1());
	P.clear();
	mult_SPD(X,s, P);
	return P;
}

inline SymMatrix mult_SPD (const RowMatrix& X)
/*
 * Symmetric Positive (Semi) Definate multiply:	return X*X'
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
 * Symmetric Positive (Semi) Definate multiply:	return X'*diag_matrix(s)*X
 */
{
	SymMatrix P(X.size2(),X.size2());
	P.clear();
	mult_SPDT(X,s, P);
	return P;
}

inline SymMatrix mult_SPDT (const RowMatrix& X)
/*
 * Symmetric Positive (Semi) Definate multiply:	return X'*X
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

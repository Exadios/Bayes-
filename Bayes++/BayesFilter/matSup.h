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
void forceSymmetric (Matrix &M, bool bUpperToLower = false);

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


}//namespace

#endif

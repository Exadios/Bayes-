#ifndef _BAYES_FILTER_MATRIX_SUPPORT
#define _BAYES_FILTER_MATRIX_SUPPORT


/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
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
 * assertion support
 */
#ifdef _DEBUG
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

RowMatrix::value_type UdUfactor_varient1 (RowMatrix& M, Subscript n);
RowMatrix::value_type UdUfactor_varient2 (RowMatrix& M, Subscript n);
inline RowMatrix::value_type UdUfactor (RowMatrix& M, Subscript n)
{	return UdUfactor_varient2(M,n);
}
LTriMatrix::value_type LdLfactor (LTriMatrix& M, Subscript n);
UTriMatrix::value_type UCfactor (UTriMatrix& M, Subscript n);

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
#ifndef BAYESFILTER_MATRIX_NOTRI
inline void Lzero (UTriMatrix& U) {}
inline void Uzero (LTriMatrix& M) {}
#endif

/*
 * Functions using UdU factorisation:
 *  inverse of Positive Definate matrix
 */
SymMatrix::value_type UdUinversePD (SymMatrix& M);
SymMatrix::value_type UdUinversePD (SymMatrix& M, SymMatrix::value_type& detM);
SymMatrix::value_type UdUinversePD (SymMatrix& MI, const SymMatrix& M);
SymMatrix::value_type UdUinversePD (SymMatrix& MI, SymMatrix::value_type& detM, const SymMatrix& M);


}//namespace

#endif

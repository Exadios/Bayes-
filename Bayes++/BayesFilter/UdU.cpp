/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Linear algebra support functions for filter classes
 * Cholesky and Modified Cholesky factorisations
 *
 * UdU' and LdL' factorisations of Positive semi-definate matrices. Where
 *  U is unit upper triangular
 *  d is diagonal
 *  L is unit lower triangular
 * Storage
 *  UD(RowMatrix) format of UdU' factor
 *   strict_upper_triangle(UD) = strict_upper_triangle(U), diagonal(UD) = d, strict_lower_triangle(UD) ignored or zeroed
 *  LD(LTriMatrix) format of LdL' factor 
 *   strict_lower_triangle(LD) = strict_lower_triangle(L), diagonal(LD) = d, strict_upper_triangle(LD) ignored or zeroed
 */
#include "bayesFlt.hpp"
#include "matSup.hpp"

#include <cassert>
#include <cmath>


/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{


template <class V>
inline typename V::value_type rcond_internal (const V& D)
/*
 * Estimate the reciprocal condition number of the diagonal vector D
 *
 * The Condition Number is defined from a matrix norm.
 *  Choose max element of D as the norm of the original matrix.
 *  Assume this norm for inverse matrix is min element D.
 * Using the D factor is fast and simple and avoids computing any squares.
 *
 * Note:
 *  Defined to be 0 for semi-definate or zero or empty matrix
 *  Defined to be <0 for negative matrix (D element a value  < 0)
 *
 *  A negative matrix may also be due to errors in the original matrix resulting in
 *   a factorisation producing special values in D (e.g. -infinity,NaN etc)
 *  By definition rcond <= 1 as min<=max
 */
{
	// Special case an empty vector
	const size_t n = D.size();
	if (n == 0)
		return 0;

	Vec::value_type rcond, mind = D[0], maxd = 0;

	for (size_t i = 0; i < n; ++i)	{
		Vec::value_type d = D[i];
		if (d < mind) mind = d;
		if (d > maxd) maxd = d;
	}
	// mind is NaN, define as negative matrix
	if (!(mind <= maxd))
		return -1;

	if (maxd == 0)	{	// avoid division by zero, for zero or negative matrix
		rcond = mind;	// mind may be zero of less  depending on D
	}
	else {
		// Rcond from min/max norm
		rcond = mind / maxd;
	}
	assert (rcond <= 1.);
	return rcond;
}

Vec::value_type rcond_vec (const Vec& D)
/*
 * Estimate the reciprocal condition number of the diagonal vector D
 * Note: diagonal(D) == UD_factor(diagonal(D)) so this rcond of both an original D and it factor
 * Equivilent to UdUrcond(diagonal(D))
 */
{
	return rcond_internal (D);
}

RowMatrix::value_type UdUrcond (const RowMatrix& UD)
/*
 * Estimate the reciprocal condition number of the original PSD
 * matrix for which UD is the factor UdU' or LdL'
 */
{
	assert (UD.size1() == UD.size2());
	return rcond_internal (diag(UD));
}

RowMatrix::value_type UdUrcond (const RowMatrix& UD, size_t n)
/*
 * As above but only first n elements are used
 */
{
	return rcond_internal (diag(UD,n));
}

UTriMatrix::value_type UCrcond (const UTriMatrix& UC)
/*
 * Estimate the reciprocal condition number of the original PSD
 * matrix for which U is the factor UU'
 *
 * See UdUrcond(RowMatrix) for definition and return value
 */
{
	assert (UC.size1() == UC.size2());
	Float rcond = rcond_internal (diag(UC));
	// Square to get rcond of original matrix, take care to propogate rcond's sign!
	if (rcond < 0)
		return -(rcond*rcond);
	else
		return rcond*rcond;
}

inline UTriMatrix::value_type UCrcond (const UTriMatrix& UC, size_t n)
/*
 * As above but for use by UCfactor functions where only first n elements are used
 */
{
	Float rcond = rcond_internal (diag(UC,n));
	// Square to get rcond of original matrix, take care to propogate rcond's sign!
	if (rcond < 0)
		return -(rcond*rcond);
	else
		return rcond*rcond;
}


RowMatrix::value_type UdUdet (const RowMatrix& UD)
/*
 * Compute the determinant of the original PSD
 * matrix for which UD is the factor UdU' or LdL'
 * Result comes directly from determinant of diagonal in triangular matrices
 *  Defined to be 1 for 0 size UD
 */
{
	const size_t n = UD.size1();
	assert (n == UD.size2());
	RowMatrix::value_type det = 1.;
	for (size_t i = 0; i < n; ++i)	{
		det *= UD(i,i);
	}
	return det;
}


RowMatrix::value_type UdUfactor_variant1 (RowMatrix& M, size_t n)
/*
 * In place Modified upper triangular Cholesky factor of a
 *  Positive definate or semi-definate matrix M
 * Reference: A+G p.218 Upper cholesky algorithm modified for UdU'
 *  Numerical stability may not be as good as M(k,i) is updated from previous results
 * Strict lower triangle of M is ignored in computation
 *
 * Input: M, n=last size_t to be included in factorisation
 * Output: M as UdU' factor
 *    strict_upper_triangle(M) = strict_upper_triangle(U)
 *    diagonal(M) = d
 *    strict_lower_triangle(M) is unmodifed
 * Return:
 *    reciprocal condition number, -1 if negative, 0 if semi-definate (including zero)
 */
{
	size_t i,j,k;
	RowMatrix::value_type e, d;

	if (n > 0)
	{
		j = n-1;
		do {
			d = M(j,j);

			// Diagonal element
			if (d > 0)
			{	// Positive definate
				d = 1 / d;

				for (i = 0; i < j; ++i)
				{
					e = M(i,j);
					M(i,j) = d*e;
					for (k = 0; k <= i; ++k)
					{
						RowMatrix::Row Mk(M,k);
						Mk[i] -= e*Mk[j];
					}
				}
			}
			else if (d == 0)
			{	// Possibly Semidefinate, check not negative
				for (i = 0; i < j; ++i)
				{
					if (M(i,j) != 0)
						goto Negative;
				}
			}
			else
			{	// Negative
				goto Negative;
			}
		} while (j-- > 0);
	}

	// Estimate the reciprocal condition number
	return rcond_internal (diag(M,n));

Negative:
   return -1.;
}


RowMatrix::value_type UdUfactor_variant2 (RowMatrix& M, size_t n)
/*
 * In place Modified upper triangular Cholesky factor of a
 *  Positive definate or semi-definate matrix M
 * Reference: A+G p.219 right side of table
 * Strict lower triangle of M is ignored in computation
 *
 * Input: M, n=last size_t to be included in factorisation
 * Output: M as UdU' factor
 *    strict_upper_triangle(M) = strict_upper_triangle(U)
 *    diagonal(M) = d
 *    strict_lower_triangle(M) is unmodifed
 * Return:
 *    reciprocal condition number, -1 if negative, 0 if semi-definate (including zero)
 */
{
	size_t i,j,k;
	RowMatrix::value_type e, d;

	if (n > 0)
	{
		j = n-1;
		do {
			RowMatrix::Row Mj(M,j);
			d = Mj[j];

			// Diagonal element
			if (d > 0)
			{	// Positive definate
				i = j;
				do
				{
					RowMatrix::Row Mi(M,i);
					e = Mi[j];
					for (k = j+1; k < n; ++k)
					{
						e -= Mi[k]*M(k,k)*Mj[k];
					}
					if (i == j) {
						Mi[j] = d = e;		// Diagonal element
					}
					else {
						Mi[j] = e / d;
					}
				} while (i-- > 0);
			}
			else if (d == 0)
			{	// Possibly Semidefinate, check not negative, whole row must be identically zero
				for (k = j+1; k < n; ++k)
				{
					if (Mj[k] != 0)
						goto Negative;
				}
			}
			else
			{	// Negative
				goto Negative;
			}
		} while (j-- > 0);
	}

	// Estimate the reciprocal condition number
	return rcond_internal (diag(M,n));

Negative:
   return -1.;
}


LTriMatrix::value_type LdLfactor (LTriMatrix& M, size_t n)
/*
 * In place Modified lower triangular Cholesky factor of a
 *  Positive definate or semi-definate matrix M
 * Reference: A+G p.218 Lower cholesky algorithm modified for LdL'
 *
 * Input: M, n=last size_t to be included in factorisation
 * Output: M as LdL' factor
 *    strict_lower_triangle(M) = strict_lower_triangle(L)
 *    diagonal(M) = d
 * Return:
 *    reciprocal condition number, -1 if negative, 0 if semi-definate (including zero)
 * ISSUE: This could change to equivilient of UdUfactor_varient2
 */
{
	size_t i,j,k;
	LTriMatrix::value_type e, d;

	for (j = 0; j < n; ++j)
	{
		d = M(j,j);

		// Diagonal element
		if (d > 0)
		{
			// Positive definate
			d = 1 / d;

			for (i = j+1; i < n; ++i)
			{
				e = M(i,j);
				M(i,j) = d*e;
				for (k = i; k < n; ++k)
				{
					LTriMatrix::Row Mk(M,k);
					Mk[i] -= e*Mk[j];
				}
			}
		}
		else if (d == 0)
		{
			// Possibly Semidefinate, check not negative
			for (i = j+1; i < n; ++i)
			{
				if (M(i,j) != 0)
					goto Negative;
			}
		}
		else
		{
			// Negative
			goto Negative;
		}
	}

	// Estimate the reciprocal condition number
	return rcond_internal (diag(M,n));

Negative:
	return -1.;
}


UTriMatrix::value_type UCfactor (UTriMatrix& M, size_t n)
/*
 * In place Upper triangular Cholesky factor of a
 *  Positive definate or semi-definate matrix M
 * Reference: A+G p.218
 * Strict lower triangle of M is ignored in computation
 *
 * Input: M, n=last size_t to be included in factorisation
 * Output: M as UC*UC' factor
 *    upper_triangle(M) = UC
 * Return:
 *    reciprocal condition number, -1 if negative, 0 if semi-definate (including zero)
 */
{
	size_t i,j,k;
	UTriMatrix::value_type e, d;

	if (n > 0)
	{
		j = n-1;
		do {
			d = M(j,j);

			// Diagonal element
			if (d > 0)
			{
				// Positive definate
				d = std::sqrt(d);
				M(j,j) = d;
				d = 1 / d;

				for (i = 0; i < j; ++i)
				{
					e = d*M(i,j);
					M(i,j) = e;
					for (k = 0; k <= i; ++k)
					{
						UTriMatrix::Row Mk(M,k);
						Mk[i] -= e*Mk[j];
					}
				}
			}
			else if (d == 0)
			{
				// Possibly Semidefinate, check not negative
				for (i = 0; i < j; ++i)
				{
					if (M(i,j) != 0)
						goto Negative;
				}
			}
			else
			{
				// Negative
				goto Negative;
			}
		} while (j-- > 0);
	}

	// Estimate the reciprocal condition number
	return UCrcond (M,n);

Negative:
   return -1;
}


RowMatrix::value_type UdUfactor (RowMatrix& UD, const SymMatrix& M)
/*
 * Modified upper triangular Cholesky factor of a
 * Positive definate or Semi-definate Matrix M
 * Wraps UdUfactor for non in place factorisation
 * Output:
 *    UD the UdU' factorisation of M with strict lower triangle zero
 * Return:
 *    see in-place UdUfactor
 */
{
	UD.assign (M);
	RowMatrix::value_type rcond = UdUfactor (UD, M.size1());

	Lzero (UD);	// Zero lower triangle ignored by UdUfactor
	return rcond;
}


LTriMatrix::value_type LdLfactor (LTriMatrix& LD, const SymMatrix& M)
/*
 * Modified lower triangular Cholesky factor of a
 * Positive definate or Semi-definate Matrix M
 * Wraps LdLfactor for non in place factorisation
 * Output:
 *    LD the LdL' factorisation of M
 * Return:
 *    see in-place LdLfactor
 */
{
	LD.assign (M);
	LTriMatrix::value_type rcond = LdLfactor (LD, M.size1());

	return rcond;
}


UTriMatrix::value_type UCfactor (UTriMatrix& UC, const SymMatrix& M)
/*
 * Upper triangular Cholesky factor of a
 * Positive definate or Semi-definate Matrix M
 * Wraps UCfactor for non in place factorisation
 * Output:
 *    UC the UC*UC' factorisation of M
 * Return:
 *    see in-place UCfactor
 */
{
	UC.assign (UTriMatrix(M));
	UTriMatrix::value_type rcond = UCfactor (UC, M.size1());

	return rcond;
}



bool UdUinverse (RowMatrix& UD)
/*
 * In-place (destructive) inversion of diagonal and unit upper triangular matrices in UD
 * BE VERY CAREFUL THIS IS NOT THE INVERSE OF UD
 *  Inversion on d and U is seperate: inv(U)*inv(d)*inv(U') = inv(U'dU) NOT EQUAL inv(UdU')
 * Lower triangle of UD is ignored and unmodified
 * Only diagonal part d can be singular (zero elements), inverse is computed of all elements other then singular
 * Reference: A+G p.223
 *
 * Output:
 *    UD: inv(U), inv(d)
 * Return:
 *    singularity (of d), true iff d has a zero element
 */
{
	size_t i,j,k;
	const size_t n = UD.size1();
	assert (n == UD.size2());

	// Invert U in place
	if (n > 1)
	{
		i = n-2;
		do {
			RowMatrix::Row UDi(UD,i);
			for (j = n-1; j > i; --j)
			{
				RowMatrix::value_type UDij = - UDi[j];
				for (k = i+1; k < j; ++k)
					UDij -= UDi[k] * UD(k,j);
				UDi[j] = UDij;
			}
		} while (i-- > 0);
	}

	// Invert d in place
	bool singular = false;
	for (i = 0; i < n; ++i)
	{
		// Detect singular element
		if (UD(i,i) != 0)
			UD(i,i) = 1 / UD(i,i);
		else
			singular = true;
	}

	return singular;
}


bool UTinverse (UTriMatrix& U)
/*
 * In-place (destructive) inversion of upper triangular matrix in U
 *
 * Output:
 *    U: inv(U)
 * Return:
 *    singularity (of U), true iff diagonal of U has a zero element
 */
{
	const size_t size = U.size1();
	assert (size == U.size2());

	bool singular = false;
	UTriMatrix UI(size,size);	// ISSUE Temporay matrix
	identity(UI);

	// ISSUE Cannot use ublas::inplace_solve (U, UI, ublas::upper_tag());
	// It only checks for singularity in debug build
	// Reimplement the general algorithm with a few tweaks

	for (UTriMatrix::difference_type n = size - 1; n >= 0; -- n) {
		UTriMatrix::value_type d = U(n, n);
		if (d == 0) {
			singular = true;
			break;
		}
		for (UTriMatrix::difference_type l = size - 1; l >= 0; -- l) {
			UTriMatrix::value_type t;
			if (n <= l)		// ISSUE uBLAS may not support indexing a 0 element
				t = UI(n, l) / d;
			else
				t = 0;
			if (t != 0)		// Optimises for zeros in UI and prevents assignments to 0 elements of packed matrices
			{
				UI(n, l) = t;
				UTriMatrix::const_reverse_iterator1 it1e1 (U.find_first1 (1, n, n));
				UTriMatrix::const_reverse_iterator1 it1e1_rend (U.find_first1 (1, 0, n));
				while (it1e1 != it1e1_rend) {
					UI(it1e1.index1 (), l) -= *it1e1 * t;
					++ it1e1;
				}
			}
		}
	}

	U = UI;
	return singular;
}


void UdUrecompose_transpose (RowMatrix& M)
/*
 * In-place recomposition of Symetric matrix from U'dU factor store in UD format
 *  Generally used for recomposing result of UdUinverse
 * Note definateness of result depends purely on diagonal(M)
 *  i.e. if d is positive definate (>0) then result is positive definate
 * Reference: A+G p.223
 * In place computation uses simple structure of solution due to triangular zero elements
 *  Defn: R = (U' d) row i , C = U column j   -> M(i,j) = R dot C;
 *  However M(i,j) only dependant R(k<=i), C(k<=j) due to zeros
 *  Therefore in place multiple sequences such k < i <= j
 * Input:
 *    M - U'dU factorisation (UD format)
 * Output:
 *    M - U'dU recomposition (symetric)
 */
{
	size_t i,j,k;
	const size_t n = M.size1();
	assert (n == M.size2());

	// Recompose M = (U'dU) in place
	if (n > 0)
	{
		i = n-1;
		do {
			RowMatrix::Row Mi(M,i);
			// (U' d) row i of lower triangle from upper trinagle
			for (j = 0; j < i; ++j)
				Mi[j] = M(j,i) * M(j,j);
			// (U' d) U in place
			j = n-1;
			do { // j>=i
				// Compute matrix product (U'd) row i * U col j
				RowMatrix::value_type Mij = Mi[j];
				if (j > i)					// Optimised handling of 1 in U
					Mij *= Mi[i];
				for (k = 0; k < i; ++k)		// Inner loop k < i <=j, only strict triangular elements
					Mij += Mi[k] * M(k,j);		// M(i,k) element of U'd, M(k,j) element of U
				M(j,i) = Mi[j] = Mij;
			} while (j-- > i);
		} while (i-- > 0);
	}
}


void UdUrecompose (RowMatrix& M)
/*
 * In-place recomposition of Symetric matrix from UdU' factor store in UD format
 *  See UdUrecompose_transpose()
 * Input:
 *    M - UdU' factorisation (UD format)
 * Output:
 *    M - UdU' recomposition (symetric)
 */
{
	size_t i,j,k;
	const size_t n = M.size1();
	assert (n == M.size2());

	// Recompose M = (UdU') in place
	for (i = 0; i < n; ++i)
	{
		RowMatrix::Row Mi(M,i);
		// (d U') col i of lower triangle from upper trinagle
		for (j = i+1; j < n; ++j) {
			RowMatrix::Row Mj(M,j);
			Mj[i] = M(i,j) * Mj[j];
		}
		// U (d U') in place
		for (j = 0; j <= i; ++j)	// j<=i
		{
			// Compute matrix product (U'd) row i * U col j
			RowMatrix::value_type Mij = Mi[j];
			if (j > i)					// Optimised handling of 1 in U
				Mij *= Mi[i];
			for (k = i+1; k < n; ++k)		// Inner loop k > i >=j, only strict triangular elements
				Mij += Mi[k] * M(k,j);		// M(i,k) element of U'd, M(k,j) element of U
			M(j,i) = Mi[j] = Mij;
		}
	}
}


void Lzero (RowMatrix& M)
/*
 * Zero strict lower triangle of Matrix
 */
{
	size_t i,j;
	const size_t n = M.size1();
	assert (n == M.size2());
	for (i = 1; i < n; ++i)
	{
		RowMatrix::Row Ui(M,i);
		for (j = 0; j < i; ++j)
		{
			Ui[j] = 0;
		}
	}
}

void Uzero (RowMatrix& M)
/*
 * Zero strict upper triangle of Matrix
 */
{
	size_t i,j;
	const size_t n = M.size1();
	assert (n == M.size2());
	for (i = 0; i < n; ++i)
	{
		RowMatrix::Row Li(M,i);
		for (j = i+1; j < n; ++j)
		{
			Li[j] = 0;
		}
	}
}


void UdUfromUCholesky (RowMatrix& U)
/*
 * Convert a normal upper triangular Cholesky factor into
 * a Modified Cholesky factor.
 * Lower triangle of UD is ignored and unmodified
 * Ignores Columns with zero diagonal element
 *  Correct for zero columns i.e. UD is Cholesky factor of a PSD Matrix
 * Note: There is no inverse to this function toCholesky as square losses the sign
 *
 * Input:
 *    U Normal Cholesky factor (Upper triangular)
 * Output:
 *    U Modified Cholesky factor (UD format)
 */
{
	size_t i,j;
	const size_t n = U.size1();
	assert (n == U.size2());
	for (j = 0; j < n; ++j)
	{
		RowMatrix::value_type sd = U(j,j);
		U(j,j) = sd*sd;
					// Devide columns by square of non zero diagonal
		if (sd != 0)
		{
			for (i = 0; i < j; ++i)
			{
				U(i,j) /= sd;
			}
		}
	}
}

void UdUseperate (RowMatrix& U, Vec& d, const RowMatrix& UD)
/*
 * Extract the seperate U and d parts of the UD factorisation
 * Output:
 *    U and d parts of UD
 */
{
	size_t i,j;
	const size_t n = UD.size1();
	assert (n == UD.size2());

	for (j = 0; j < n; ++j)
	{
					// Extract d and set diagonal to 1
		d[j] = UD(j,j);
		RowMatrix::Row Uj(U,j);
		Uj[j] = 1.;
		for (i = 0; i < j; ++i)
		{
			U(i,j) = UD(i,j);
			// Zero lower triangle of U
			Uj[i] = 0;
		}
	}
}

/*
 * Function built using UdU factorisation
 */

SymMatrix::value_type UdUinversePD (SymMatrix& M)
/*
 * inverse of Positive Definate matrix
 * Input:
 *     M is a symetric matrix
 * Output:
 *     M inverse of M, only updated if return value >0
 * Return:
 *     reciprocal condition number, -1 if negative, 0 if semi-definate (including zero)
 */
{
					// Abuse as a RowMatrix
	RowMatrix& M_Matrix = M.asRowMatrix();
	SymMatrix::value_type rcond = UdUfactor (M_Matrix, M.size1());
	// Only invert and recompose if PD
	if (rcond > 0) {
		(void)UdUinverse (M_Matrix);
		UdUrecompose_transpose (M_Matrix);
	}
	return rcond;
}

SymMatrix::value_type UdUinversePD (SymMatrix& M, SymMatrix::value_type& detM)
/*
 * As above but also computes determinant of original M if M is PSD
 */
{
					// Abuse as a RowMatrix
	RowMatrix& M_Matrix = M.asRowMatrix();
	SymMatrix::value_type rcond = UdUfactor (M_Matrix, M.size1());
	// Only invert and recompose if PD
	if (rcond > 0) {
		detM = UdUdet(M_Matrix);
		(void)UdUinverse (M_Matrix);
		UdUrecompose_transpose (M_Matrix);
	}
	return rcond;
}


SymMatrix::value_type UdUinversePD (SymMatrix& MI, const SymMatrix& M)
/*
 * inverse of Positive Definate matrix
 * Input:
 *    M is a symetric matrix
 * Output:
 *    MI inverse of M, only valid if return value >0
 * Return:
 *    reciprocal condition number, -1 if negative, 0 if semi-definate (including zero)
 */
{
	MI = M;
					// Abuse as a RowMatrix
	RowMatrix& MI_Matrix = MI.asRowMatrix();
	SymMatrix::value_type rcond = UdUfactor (MI_Matrix, MI.size1());
	// Only invert and recompose if PD
	if (rcond > 0) {
		(void)UdUinverse (MI_Matrix);
		UdUrecompose_transpose (MI_Matrix);
	}
	return rcond;
}

SymMatrix::value_type UdUinversePD (SymMatrix& MI, SymMatrix::value_type& detM, const SymMatrix& M)
/*
 * As above but also computes determinant of original M if M is PSD
 */
{
	MI = M;
					// Abuse as a RowMatrix
	const RowMatrix& M_Matrix = M.asRowMatrix();
	RowMatrix& MI_Matrix = MI.asRowMatrix();
	SymMatrix::value_type rcond = UdUfactor (MI_Matrix, MI.size1());
	// Only invert and recompose if PD
	if (rcond > 0) {
		detM = UdUdet(M_Matrix);
		(void)UdUinverse (MI_Matrix);
		UdUrecompose_transpose (MI_Matrix);
	}
	return rcond;
}

}//namespace

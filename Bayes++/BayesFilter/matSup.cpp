/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */

/*
 * Matrix support functions for filter classes
 * Relies on Bayesian_filter::Bayes_filter for exception
 * thrown by internal matrix checks
 */
#include "matSup.h"
#include <cassert>

#include "bayesFlt.h"		// Require Excecptions

namespace {

template <class scalar>
inline scalar sqr(scalar x)
// Square 
{
	return x*x;
}
};//namespace


/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{
	using Bayesian_filter::Bayes_filter_exception;


#ifdef _DEBUG
void assert_isSymmetric (const SymMatrix &M)
/*
 * Assert a Matrix is Square and Symmetric
 *  Requires 'cerr' and 'assert' to about execution
 */
{
	bool bSym = isSymmetric(M);
	if (!bSym)
	{
		// Display Asymmetry
		std::cerr.flags(std::ios::scientific); std::cerr.precision(17);
		std::cerr << M;
	}
	assert(bSym);
}

void assert_isPSD (const SymMatrix &M)
/*
 * Assert a Matrix is Positive Semi Definate via the alogrithm in isPSD
 *  Requires 'cerr' and 'assert' to about execution
 */
{
	bool bPSD = isPSD(M);
	if (!bPSD) {
		// Display Non PSD
		std::cerr.flags(std::ios::scientific); std::cerr.precision(17);
		std::cerr << M;
	}
	assert(bPSD);
}
#endif


bool isPSD (const SymMatrix &M)
/*
 * Check a Matrix is both Symertic and Positive Semi Definate
 *  Creates a temporary copy
 *  
 * Numerics of Algorithm:
 *  Use UdUfactor and checks reciprocal condition number >=0
 *  Algorithms using UdUfactor will will therefore succeed if isPSD(M) is true
 *  Do NOT assume because isPSD is true that M will appear to be PSD to other numerical algorithms
 * Return:
 *  true iff M isSymetric and is PSD by the above algorithm
 */
{
	// Precondition Square and Symmetric
	if (!isSymmetric (M))
		return false;

	RowMatrix UD(M.size1(),M.size1());

	RowMatrix::value_type rcond = UdUfactor(UD, M);
	return rcond >= 0.;
}


bool isSymmetric (const SymMatrix &M)
/*
 * Check a Symmetric Matrix really is Square and Symmetric
 * Implictly also checks matrix is without IEC 559 NaN values as they are always !=
 * Return:
 *  true iff M is Square and Symmetric and without NaN 
 */
{
	// Check Square
	if (M.size1() != M.size2() ) {
		throw Bayes_filter_exception ("SymMatrix is not square");
	}

	// Check Symmetric
	bool bSym = true;
	Subscript size = M.size1();
	for (Subscript r = 0; r < size; ++r) {
		for (Subscript c = 0; c <= r; ++c) {
			if( M(r,c) != M(c,r) ) {
				bSym = false;
			}
		}
	}
	return bSym;
}


void forceSymmetric (SymMatrix &M, bool bUpperToLower)
/*
 * Force Matrix Symmetry
 *	Normally Copies ower triangle to upper or
 *   Upper to Lower triangle is specifed
 */
{
	// Check Square
	if (M.size1() != M.size2() ) {
		throw Bayes_filter_exception ("SymMatrix is not square");
	}

	Subscript size = M.size1();

	if (bUpperToLower)
	{
		// Copy Lower to Upper
		for (Subscript r = 1; r < size; ++r) {
			for (Subscript c = 0; c < r; ++c) {
				M(c,r) = M(r,c);
			}
		}
	}
	else
	{
		// Copy Upper to Lower
		for (Subscript r = 1; r < size; ++r) {
			for (Subscript c = 0; c < r; ++c) {
				M(r,c) = M(c,r);
			}
		}
	}
}



}//namespace

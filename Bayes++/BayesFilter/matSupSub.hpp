/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Matrix types for filter classes
 *  Replace this header to replace matrix support functions
 *
 * Everything in namespace Bayes_filter_matrix is intended to support the matrix storage
 * and algebra requirements of the library. Therefore the interfaces and implementation is
 * not intended to be stable.
 */

/*
 * Use the Boost uBLAS Basic Linear Algebra library
 * That is boost::numeric::ublas
 *  Thanks to  Joerg Walter, Mathias Koch for an excellent C++ library
 *
 * Provides predefined types Vec and a variety of Matrix types
 *
 * Sparse support: The macros BAYES_FILTER_(SPARSE/COMPRESSED/COORDINATE) control experimental sparse matrix support
 * These simply replace the default storage types with their sparse equivilents
 */

#include <boost/version.hpp>
#if !(BOOST_VERSION >= 103000)
#error Requires Boost 1.30.0 or later
#endif

// Element proxies have a colourful history! They are not required by Bayes++
// As of Boost 1.30.0 they do not allow mixed assignment of elements, and must be disabled
#define BOOST_UBLAS_NO_ELEMENT_PROXIES

// Required to allow all members of trinangular matrices to be accessed
#define BOOST_UBLAS_REFERENCE_CONST_MEMBER

// Singularity is a numerical issue. Thanks to Joerg we can enable runtime exceptions
#define BOOST_UBLAS_USE_EXCEPTIONS
#define BOOST_UBLAS_SINGULAR_CHECK

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#if defined(BAYES_FILTER_SPARSE) || defined(BAYES_FILTER_COMPRESSED) || defined(BAYES_FILTER_COORDINATE)
#include <map>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#define BAYES_FILTER_GAPPY
#endif


/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{
						// Allow use of ublas
namespace ublas = boost::numeric::ublas;

/*
 * Declare the value used for ALL linear algebra operations
 * Also required as the matrix/vector container value_type
 */
typedef double Float;

/*
 * uBlas base types - these will be wrapper to provide the actual vector and matrix types
 *  Symmetric types don't appear. They are defined later by adapting these base types
 */
namespace detail {
							// Dense types
typedef ublas::vector<Float> BaseDenseVector;
typedef ublas::matrix<Float, ublas::row_major> BaseDenseRowMatrix;
typedef ublas::matrix<Float, ublas::column_major> BaseDenseColMatrix;
typedef ublas::triangular_matrix<Float, ublas::upper, ublas::row_major> BaseDenseUpperTriMatrix;
typedef ublas::triangular_matrix<Float, ublas::lower, ublas::row_major> BaseDenseLowerTriMatrix;
typedef ublas::banded_matrix<Float> BaseDenseDiagMatrix;
							// Sparse types
#if defined(BAYES_FILTER_SPARSE)
typedef ublas::sparse_vector<Float, std::map<size_t,Float> > BaseSparseVector;
typedef ublas::sparse_matrix<Float, ublas::row_major, std::map<size_t,Float> > BaseSparseRowMatrix;
typedef ublas::sparse_matrix<Float, ublas::column_major, std::map<size_t,Float> > BaseSparseColMatrix;
							// OR Compressed types
#elif defined(BAYES_FILTER_COMPRESSED)
typedef ublas::compressed_vector<Float> BaseSparseVector;
typedef ublas::compressed_matrix<Float, ublas::row_major> BaseSparseRowMatrix;
typedef ublas::compressed_matrix<Float, ublas::column_major> BaseSparseColMatrix;
							// OR Coordinate types
#elif defined(BAYES_FILTER_COORDINATE)
typedef ublas::coordinate_vector<Float> BaseSparseVector;
typedef ublas::coordinate_matrix<Float, ublas::row_major> BaseSparseRowMatrix;
typedef ublas::coordinate_matrix<Float, ublas::column_major> BaseSparseColMatrix;
#endif

							// Default types Dense or Gappy
#ifndef BAYES_FILTER_GAPPY
typedef BaseDenseVector BaseVector;
typedef BaseDenseRowMatrix BaseRowMatrix;
typedef BaseDenseColMatrix BaseColMatrix;
typedef BaseDenseUpperTriMatrix BaseUpperTriMatrix;
typedef BaseDenseLowerTriMatrix BaseLowerTriMatrix;
typedef BaseDenseDiagMatrix BaseDiagMatrix;
#else
typedef BaseSparseVector BaseVector;
typedef BaseSparseRowMatrix BaseRowMatrix;
typedef BaseSparseColMatrix BaseColMatrix;
typedef BaseDenseUpperTriMatrix BaseUpperTriMatrix;		// No sparse triangular or banded
typedef BaseDenseLowerTriMatrix BaseLowerTriMatrix;
typedef BaseDenseDiagMatrix BaseDiagMatrix;
#endif

}

}//namespace

/*
 * Common type independant uBlas interface
 */
#include "uBLASmatrix.hpp"

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Matrix support functions sub header for filter classes
 * Replace this header to replace matrix support functions
 */

/*
 * Matrix types for filter classes
 * Implemented using boost::numeric::ublas uBblas Basic Linear Algebra library
 *  Provides predefined types Vec and a variety of Matrix types
 *
 * Everything in namespace Bayes_filter_matrix is intended to support the matrix storage
 * and algebra requirements of the library. Therefore the interfaces and implementation is
 * not intended to be stable.
 *
 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

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
 *  We require static type conversion between RowMatrix and SymMatrix
 *  This requires they both use the same dense represenation. Therefore
 *  we use a symmetric_adaptor to provide the base for symmetric matrices.
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
typedef ublas::sparse_vector<Float> BaseSparseVector;
typedef ublas::sparse_matrix<Float, ublas::row_major> BaseSparseRowMatrix;
typedef ublas::sparse_matrix<Float, ublas::column_major> BaseSparseColMatrix;

							// Default types
typedef BaseSparseVector BaseVector;
typedef BaseSparseRowMatrix BaseRowMatrix;
typedef BaseSparseColMatrix BaseColMatrix;
typedef BaseSparseRowMatrix BaseUpperTriMatrix;
typedef BaseSparseRowMatrix BaseLowerTriMatrix;
typedef BaseSparseRowMatrix BaseDiagMatrix;
/*
typedef BaseDenseVector BaseVector;
typedef BaseDenseRowMatrix BaseRowMatrix;
typedef BaseDenseColMatrix BaseColMatrix;
typedef BaseDenseRowMatrix BaseUpperTriMatrix;
typedef BaseDenseRowMatrix BaseLowerTriMatrix;
typedef BaseDenseRowMatrix BaseDiagMatrix;*/

}

}//namespace

/*
 * Common type independant uBlas interface
 */
#include "uBLASmatrix.hpp"

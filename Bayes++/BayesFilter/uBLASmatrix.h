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
 * Implemented using boost::numeric::ublas uBLAS Basic Linear Algebra library
 *  Provides predefined types Vec and a variety of Matrix types
 *
 * Everything in namespace Bayes_filter_matrix is intended to support the matrix storage
 * and algebra requirements of the library. Therefore the interfaces and implementation is
 * not intended to be stable. Nor is this a general purpose adapator for uBLAS
 *
 * Note on row_major matrices
 *  The implementation uses row major extensively. The computation of symetric products P*P' is 
 *  most efficient with row operations. These products are used extensively so the default
 *  is to use row_major matrices
 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>

/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{
						// Allow use of ublas and a few functions in own namespace
namespace ublas = boost::numeric::ublas;
using ublas::row;
using ublas::column;
using ublas::trans;
using ublas::prod;			// These do not applie to the templated prod<temp> funtions
using ublas::inner_prod;
using ublas::outer_prod;

enum EmptyTag {Empty};	// Tag type used for empty matrix constructor



/*
 * Filter Vec type
 */
class Vec: public ublas::vector<double>
{
	typedef ublas::vector<double> VecBase;
public:
	// No Default Constructor. Empty creation is very error prone
	explicit Vec(EmptyTag) : VecBase()
	{}	// Empty constructor
	explicit Vec(size_t s) : VecBase(s)
	{}	// Normal sized constructor
	template <class E>
	Vec(const ublas::vector_expression<E>& e) : VecBase(e)
	{}	// vector_expression conversion constructor

	Vec& operator=(const Vec& r)
	{	// Vector assignment; independant
		assign(r);
		return *this;
	}
	template <class Base>
	Vec& operator=(const ublas::matrix_row<Base>& r)
	{	// Matrix row assignment; independant
		assign(r);
		return *this;
	}
	template <class Base>
	Vec& operator=(const ublas::matrix_column<Base>& r)
	{	// Matrix column assignment; independant
		assign(r);
		return *this;
	}
	template <class E>
	Vec& operator=(const ublas::vector_expression<E>& r)
	{	// Expression assignment, may be dependant on r
		VecBase::operator=(r);
		return *this;
	}

	const ublas::vector_range<VecBase> operator()(size_t b, size_t e) const
	{	// Range selection operator
		return ublas::vector_range<VecBase>(*const_cast<Vec*>(this), ublas::range(b,e));
	}
	ublas::vector_range<VecBase> operator()(size_t b, size_t e)
	{	// Range selection operator
		return ublas::vector_range<VecBase>(*this, ublas::range(b,e));
	}
};


namespace detail {		// Lots of implementation detail

/*
 * Filter Matrix class template. Augmentation for uBlas MatrixBase
 */
template <class MatrixBase>
class FMMatrix : public MatrixBase
{
public:
	typedef typename MatrixBase::value_type value_type;

	// No Default Constructor. Empty creation is very error prone
	explicit FMMatrix(EmptyTag) : MatrixBase()
	{}	// Empty constructor
	FMMatrix(size_t r, size_t c) : MatrixBase(r,c)
	{}	// Normal sized constructor
	FMMatrix(const FMMatrix& c) : MatrixBase(static_cast<const MatrixBase&>(c))
	{}	// Copy constructor
	template <class E>
	FMMatrix(const ublas::matrix_expression<E>& e) : MatrixBase(e)
	{}	// matrix_expression conversion constructor

	FMMatrix& operator=(const FMMatrix& r)
	{	// Matrix assignment
		assign(r);
		return *this;
	}
	template <class E>
	FMMatrix& operator=(const ublas::matrix_expression<E>& r)
	{	// Expression assignment, may be dependant on r
		MatrixBase::operator=(r);
		return *this;
	}

	// Row,Column vector proxies
	typedef ublas::matrix_row<FMMatrix> Row;
	typedef const ublas::matrix_row<const FMMatrix> const_Row;
	typedef ublas::matrix_column<FMMatrix> Column;
	typedef const ublas::matrix_column<const FMMatrix> const_Column;

	// Vector proxies from iterators - static members dependant on MatrixBase type
	// ri() returns container associated with iterator. static_cast required as typeof(ri()) may not be MM
	static ublas::matrix_row<MatrixBase> rowi(const typename MatrixBase::iterator1& ri)
	{
		typedef MatrixBase MM;
		return ublas::matrix_row<MM>(static_cast<MM&>(ri()), ri.index1());
	}
	static ublas::matrix_row<const MatrixBase> rowi(const typename MatrixBase::const_iterator1& ri)
	{
		typedef const MatrixBase MM;
		return ublas::matrix_row<MM>(static_cast<MM&>(ri()), ri.index1());
	}
	static ublas::matrix_column<MatrixBase> columni(const typename MatrixBase::iterator2& ci)
	{
		typedef MatrixBase MM;
		return ublas::matrix_column<MM>(static_cast<MM&>(ci()), ci.index2());
	}
	static ublas::matrix_column<const MatrixBase> columni(const typename MatrixBase::const_iterator2& ci)
	{
		typedef const MatrixBase MM;
		return ublas::matrix_column<MM>(static_cast<MM&>(ci()), ci.index2());
	}

	// Sub matrix/vector helpers
	ublas::matrix_range<const MatrixBase>
	sub_matrix(size_t s1, size_t e1, size_t s2, size_t e2) const
	{
		using namespace ublas;
		return matrix_range<const MatrixBase> (*this, range(s1,e1), range(s2,e2));
	}
	ublas::matrix_range<MatrixBase>
	sub_matrix(size_t s1, size_t e1, size_t s2, size_t e2)
	{
		using namespace ublas;
		return matrix_range<MatrixBase> (*this, range(s1,e1), range(s2,e2));
	}

		// ISSUE uBLAS currently has bugs when matrix_vector_slice is require to represent a sub column
#ifndef BAYESFILTER_UBLAS_SLICE_OK
		// Workaround using a vector_range of a matrix_column
private:
	typedef ublas::vector_range<ublas::matrix_column<MatrixBase> > special_matrix_vector_slice_base;
	struct special_matrix_vector_slice : special_matrix_vector_slice_base {
		special_matrix_vector_slice(MatrixBase& mb, size_t s1, size_t e1, size_t s2) :
			special_matrix_vector_slice_base(col,ublas::range(s1,e1)), col(mb, s2)
		{}
		ublas::matrix_column<MatrixBase> col;
	};
public:
	special_matrix_vector_slice
	sub_column(size_t s1, size_t e1, size_t s2)
	// Column vector s2 with rows [s1,e1)
	{
		return special_matrix_vector_slice(*this, s1,e1, s2);
	}
#else
		// For this to work requires patched uBLAS for Bayes++
	ublas::matrix_vector_slice<const MatrixBase>
	sub_column(size_t s1, size_t e1, size_t s2) const 
	// Column vector s2 with rows [s1,e1)
	{
		using namespace ublas;
		return matrix_vector_slice<const MatrixBase> (*this, slice(s1,1,e1-s1), slice(s2,0,e1-s1));
	}
	ublas::matrix_vector_slice<MatrixBase>
	sub_column(size_t s1, size_t e1, size_t s2)
	// Column vector s2 with rows [s1,e1)
	{
		using namespace ublas;
		return matrix_vector_slice<MatrixBase> (*this, slice(s1,1,e1-s1), slice(s2,0,e1-s1));
	}
#endif
};


/*
 * uBLAS Base Types
 *  We require static type converion between RowMatrix and SymMatrix
 *  This requires they both use the same dense represenation. Therefore
 *  we use a symmetric_adaptor to provide the base for symmetric matrices.
 */
typedef ublas::matrix<double, ublas::row_major> BaseRowMatrix;
typedef ublas::matrix<double, ublas::column_major> BaseColMatrix;
typedef ublas::triangular_matrix<double, ublas::upper, ublas::row_major> BaseUpperTriMatrix;
typedef ublas::triangular_matrix<double, ublas::lower, ublas::row_major> BaseLowerTriMatrix;
typedef ublas::banded_matrix<double> BaseDiagMatrix;

// Helper class for _BaseSymMatrix allow construction of BaseRowMatrix (rm) before symmertic_adaptor
class BaseSymMatrix;
class RMConstruct
{
	BaseRowMatrix rm;
	friend class BaseSymMatrix;
	RMConstruct () : rm()
	{}
	RMConstruct (BaseRowMatrix::size_type size1, BaseRowMatrix::size_type size2) : rm(size1,size2)
	{}
	RMConstruct (const BaseRowMatrix& r) : rm(r)
	{}
	template <class E>
	RMConstruct (const ublas::matrix_expression<E>& e) : rm(e)
	{}
};

// Symmetric matrix base using addapted row matrix base
class BaseSymMatrix : private RMConstruct, public ublas::symmetric_adaptor<BaseRowMatrix, ublas::upper>
{
	typedef ublas::symmetric_adaptor<BaseRowMatrix, ublas::upper> SymAdap;
public:
	BaseSymMatrix () : RMConstruct(), SymAdap(rm)
	{}
	BaseSymMatrix (size_type nsize1, size_type nsize2) : RMConstruct(nsize1,nsize2), SymAdap(rm)
	{}
	explicit BaseSymMatrix (const BaseSymMatrix& r) : RMConstruct(reinterpret_cast<const BaseRowMatrix&>(r)), SymAdap(rm)
	{}
	// Explict construction referencing a BaseRowMatrix
	template <class E>
	explicit BaseSymMatrix (const ublas::matrix_expression<E>& e) : RMConstruct(e), SymAdap(rm)
	{}
	// Explict matrix_expression conversion constructor

	template <class E>
	BaseSymMatrix& operator=(const ublas::matrix_expression<E>& r)
	{
		SymAdap::operator=(r);
		return *this;
	}

	// Conversions straight to a FMMatrix, equivilent to a RowMatrix
	const FMMatrix<BaseRowMatrix>& asRowMatrix() const
	{
		return static_cast<const FMMatrix<BaseRowMatrix>& >(rm);
	}
	FMMatrix<BaseRowMatrix>& asRowMatrix()
	{	
		return static_cast<FMMatrix<BaseRowMatrix>& >(rm);
	}

	// Matrix storage members
	void clear()
	{	rm.clear();
	}
	void resize(size_type nsize1, size_type nsize2)
	{
		rm.resize(nsize1, nsize2);
	}
};

}//namespace detail

/*
 * Define Filter Matrix types
 *  Finally the definitions !
 */
using detail::FMMatrix;	// Matrix template class for template parameter matching
typedef FMMatrix<detail::BaseRowMatrix> RowMatrix;
typedef RowMatrix Matrix;
typedef FMMatrix<detail::BaseColMatrix> ColMatrix;
typedef FMMatrix<detail::BaseSymMatrix> SymMatrix;
typedef FMMatrix<detail::BaseUpperTriMatrix> UTriMatrix;
typedef FMMatrix<detail::BaseLowerTriMatrix> LTriMatrix;
typedef FMMatrix<detail::BaseDiagMatrix> DiagMatrix;

/*
 * Type conversions helper functions
 */
inline RowMatrix& asRowMatrix(SymMatrix& M)
{
	return M.asRowMatrix();
}
inline const RowMatrix& asRowMatrix(const SymMatrix& M)
{
	return M.asRowMatrix();
}


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
 * Symmetric Positive (Semi) Definate multiplication: X*S*X' and X'*S*X
 *  The name is slightly misleading. The result is actually only PD if S is PD
 *  Algorithms are intended to exploit the symmerty of the result
 *  and also where possible the row by row multiplication inherent in X*X'
 * Algorithms vary in a
 * mult_SPD: Use a argument P += multiplication result
 * prod_SPD: Return an uBLAS expression template representing the product
 */

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
 * uBLAS Expression templates for X*X' and X*X'
 * There is no point in including composite X*S*X' form as they require tempories
 * Functions are only defined for the type for which the operation is efficient
 * ISSUE Although numerically symmetric there is no way the ET produced can define this property
 *  The result must be assigned to a symmetric container to exploit the symmetry
 */

template <class X>
struct prod_SPD_matrix_traits
{	// Provide ET result type for prod(X,trans(X))
	// ISSUE Although numerically symmetric there is no way this can be added to the ET produced here
	typedef ublas::matrix_unary2_traits<X, ublas::scalar_identity<typename X::value_type> >::result_type  XT_type;
	typedef ublas::matrix_matrix_binary_traits<typename X::value_type, X,
                                        XT_type::value_type, XT_type>::result_type  XXT_type;
};

template <class X>
struct prod_SPDT_matrix_traits
{	// Provide ET result type for prod(trans(X),X)
	// ISSUE Although numerically symmetric there is no way this can be added to the ET produced here
	typedef ublas::matrix_unary2_traits<X, ublas::scalar_identity<typename X::value_type> >::result_type  XT_type;
	typedef ublas::matrix_matrix_binary_traits<XT_type::value_type, XT_type,
                                        typename X::value_type, X>::result_type  XTX_type;
};

inline
prod_SPD_matrix_traits<RowMatrix>::XXT_type 
 prod_SPD (const RowMatrix& X, const SymMatrix& S, RowMatrix& XStemp)
/*
 * Symmetric Positive (Semi) Definate product: X*(X*S)'
 */
{
	return prod( X, trans(XStemp.assign(prod(X,S))) );
}

inline
prod_SPD_matrix_traits<RowMatrix>::XXT_type 
 prod_SPD (const RowMatrix& X)
/*
 * Symmetric Positive (Semi) Definate product: X*X'
 */
{
	return prod( X, trans(X) );
}

inline SymMatrix prod_SPD (const RowMatrix& X, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate product: X*diag_matrix(s)*X'
 * TODO: Define using ublas ET
 */
{
	SymMatrix P(X.size1(),X.size1());
	P.clear();
	mult_SPD(X,s, P);
	return P;
}

inline Vec::value_type prod_SPD (const Vec& x, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate product: x*diag_matrix(s)*x'
 * Optimised to exploit the symmetry of the computation of x * x' and result p
 * ISSUE: not sparse efficient
 * TODO: Define using ublas ET
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


inline
prod_SPDT_matrix_traits<ColMatrix>::XTX_type 
 prod_SPDT (const ColMatrix& X)
/*
 * Symmetric Positive (Semi) Definate product: X'*X
 */
{
	return prod( trans(X), X );
}

inline SymMatrix prod_SPDT (const ColMatrix& X, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate product: X'*diag_matrix(s)*X
 * TODO: Define using ublas ET
 */
{
	SymMatrix P(X.size2(),X.size2());
	P.clear();
	mult_SPDT(X,s, P);
	return P;
}

inline Vec::value_type prod_SPDT (const Vec& x, const SymMatrix& S)
/*
 * Symmetric Positive (Semi) Definate multiply:	p = x'*S*x
 */
{
	return inner_prod (x, ublas::prod(S,x) );
}


}//namespace

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Common type independant uBlas interface
 *  Should be include after base types have been defined
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

/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{
						// Allow use a few functions in own namespace
namespace ublas = boost::numeric::ublas;
using ublas::row;
using ublas::column;
using ublas::trans;
using ublas::prod;		// These do not apply to the templated prod<temp> funtions
using ublas::inner_prod;
using ublas::outer_prod;

enum EmptyTag {Empty};	// Tag type used for empty matrix constructor


namespace detail {		// Lots of implementation detail

/*
 * Filter Vec type
 */
template <class VecBase>
class FMVec : public VecBase
{
public:
	// No Default Constructor. Empty creation is very error prone
	explicit FMVec(EmptyTag) : VecBase()
	{}	// Empty constructor
	explicit FMVec(size_t s) : VecBase(s)
	{}	// Normal sized constructor
	template <class E>
	FMVec(const ublas::vector_expression<E>& e) : VecBase(e)
	{}	// vector_expression conversion constructor

	FMVec& operator=(const FMVec& r)
	{	// Vector assignment; independant
		assign(r);
		return *this;
	}
	template <class Base>
	FMVec& operator=(const ublas::matrix_row<Base>& r)
	{	// Matrix row assignment; independant
assert(false);
		return *this;
	}
	template <class Base>
	FMVec& operator=(const ublas::matrix_column<Base>& r)
	{	// Matrix column assignment; independant
assert(false);
		return *this;
	}
	template <class E>
	FMVec& operator=(const ublas::vector_expression<E>& r)
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
	typedef ublas::matrix_row<MatrixBase> Row;
	typedef const ublas::matrix_row<const MatrixBase> const_Row;
	typedef ublas::matrix_column<MatrixBase> Column;
	typedef const ublas::matrix_column<const MatrixBase> const_Column;

	// Vector proxies from iterators - static members dependant on MatrixBase type
	// ri() returns container associated with iterator. static_cast required as typeof(ri()) may not be typeof(MM)
	static Row rowi(const typename MatrixBase::iterator1& ri)
	{
		return Row(static_cast<MatrixBase&>(ri()), ri.index1());
	}
	static const_Row rowi(const typename MatrixBase::const_iterator1& ri)
	{
		return const_Row(static_cast<const MatrixBase&>(ri()), ri.index1());
	}
	static Column columni(const typename MatrixBase::iterator2& ci)
	{
		return Column(static_cast<MatrixBase&>(ci()), ci.index2());
	}
	static const_Column columni(const typename MatrixBase::const_iterator2& ci)
	{
		return const_Column(static_cast<const MatrixBase&>(ci()), ci.index2());
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


// Helper class for SymMatrixAdaptor, allow construction of BaseMatrix (rm) before symmertic_adaptor
template <class MatrixBase>
class SymMatrixAdaptor;

template <class MatrixBase>
class RMConstruct
{
	MatrixBase rm;
	friend class SymMatrixAdaptor;
	RMConstruct () : rm()
	{}
	RMConstruct (MatrixBase::size_type size1, MatrixBase::size_type size2) : rm(size1,size2)
	{}
	RMConstruct (const MatrixBase& r) : rm(r)
	{}
	template <class E>
	RMConstruct (const ublas::matrix_expression<E>& e) : rm(e)
	{}
};

// Symmetric matrix base using addaptor
template <class MatrixBase>
class SymMatrixAdaptor : private RMConstruct<MatrixBase>, public ublas::symmetric_adaptor<MatrixBase, ublas::upper>
{
	typedef ublas::symmetric_adaptor<MatrixBase, ublas::upper> SymAdap;
public:
	SymMatrixAdaptor () : RMConstruct<MatrixBase>(), SymAdap(rm)
	{}
	SymMatrixAdaptor (size_type nsize1, size_type nsize2) : RMConstruct<MatrixBase>(nsize1,nsize2), SymAdap(rm)
	{}
	explicit SymMatrixAdaptor (const SymMatrixAdaptor& r) : RMConstruct<MatrixBase>(reinterpret_cast<const MatrixBase&>(r)), SymAdap(rm)
	{}
	// Explict copy construction referencing the copy reinterpreted as a MatrixBase
	template <class E>
	explicit SymMatrixAdaptor (const ublas::matrix_expression<E>& e) : RMConstruct<MatrixBase>(e), SymAdap(rm)
	{}
	// Explict matrix_expression conversion constructor

	template <class E>
	SymMatrixAdaptor& operator=(const ublas::matrix_expression<E>& r)
	{
		SymAdap::operator=(r);
		return *this;
	}

	// Conversions straight to a FMMatrix, equivilent to a RowMatrix
	const FMMatrix<MatrixBase>& asRowMatrix() const
	{
		return static_cast<const FMMatrix<MatrixBase>& >(rm);
	}
	FMMatrix<MatrixBase>& asRowMatrix()
	{	
		return static_cast<FMMatrix<MatrixBase>& >(rm);
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
 * Define Filter Vector / Matrix types
 *  Finally the definitions !
 */
using detail::FMVec;		// Template class for template parameter matching
using detail::FMMatrix;

							// Default types
typedef FMVec<detail::BaseVector> Vec;
typedef FMMatrix<detail::BaseRowMatrix> RowMatrix;
typedef RowMatrix Matrix;
typedef FMMatrix<detail::BaseColMatrix> ColMatrix;
typedef FMMatrix<detail::SymMatrixAdaptor<detail::BaseRowMatrix> > SymMatrix;
typedef FMMatrix<detail::BaseUpperTriMatrix> UTriMatrix;
typedef FMMatrix<detail::BaseLowerTriMatrix> LTriMatrix;
typedef FMMatrix<detail::BaseDiagMatrix> DiagMatrix;

							// Explicitly dense types
typedef FMVec<detail::BaseDenseVector> DenseVec;
typedef FMMatrix<detail::BaseDenseRowMatrix> DenseRowMatrix;
typedef DenseRowMatrix DenseMatrix;
typedef FMMatrix<detail::BaseDenseColMatrix> DenseColMatrix;
typedef FMMatrix<detail::SymMatrixAdaptor<detail::BaseDenseRowMatrix> > DenseSymMatrix;
typedef FMMatrix<detail::BaseDenseUpperTriMatrix> DenseUTriMatrix;
typedef FMMatrix<detail::BaseDenseLowerTriMatrix> DenseLTriMatrix;
typedef FMMatrix<detail::BaseDenseDiagMatrix> DenseDiagMatrix;

							// Explicitly sparse types
typedef FMVec<detail::BaseSparseVector> SparseVec;
typedef FMMatrix<detail::BaseDenseRowMatrix> SparseRowMatrix;
typedef SparseRowMatrix SparseMatrix;
typedef FMMatrix<detail::BaseSparseColMatrix> SparseColMatrix;
//typedef FMMatrix<detail::SymMatrixAdaptor<detail::BaseSparseRowMatrix> > SparseSymMatrix;

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
// Clear and fill square submatrix of I as identity matrix
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
	typename MatrixX::const_iterator1 Xa = X.begin1();
	const typename MatrixX::const_iterator1 Xend = X.end1();
	typename MatrixX::const_iterator1 Xb;

	// P(a,b) = X.row(a) * X.row(b)
	for (; Xa != Xend; ++Xa)				// Iterate Rows
	{		
		MatrixX::const_Row Xav = MatrixX::rowi(Xa);
		Xb = Xa;							// Start at the row Xa only one triangle of symetric result required
		for (; Xb != Xend; ++Xb)
		{
			SymMatrix::value_type p = 0;	// Tripple vector inner product
			MatrixX::const_Row Xbv = MatrixX::rowi(Xb);
			for (si = s.begin(); si != send; ++si) {	// TODO: use iterator on Xav and Xbv
				Vec::size_type i = si.index();
				p += Xav[i] * (*si) * Xbv[i];
			}
			P(Xa.index1(),Xb.index1()) += p;
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
	typename MatrixX::const_iterator1 Xa = X.begin1();
	const typename MatrixX::const_iterator1 Xend = X.end1();
	typename MatrixX::const_iterator1 Xb;

	// P(a,b) = X.row(a) * X.row(b)
	for (; Xa != Xend; ++Xa,++Pa)			// Iterate vectors
	{		
		MatrixX::const_Row Xav = MatrixX::rowi(Xa);
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		for (; Xb != Xend; ++Xb)
		{									// Simply multiple row Xa by row Xb
			SymMatrix::value_type p = inner_prod( Xav, MatrixX::rowi(Xb) );
			P(Xa.index1(),Xb.index1()) += p;
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
	typename MatrixX::const_iterator1 Xa = X.begin1();
	const typename MatrixX::const_iterator1 Xend = X.end1();
	typename MatrixX::const_iterator1 Xb;

	// P(a,b) = X.row(a) * S * X'.col(b) = X.row(a) * S * X.row(b)
	//        = (S * X.row(a)')' * X.row(b)

	for (; Xa != Xend; ++Xa)				// Iterate Rows
	{
		stemp.assign( ublas::prod(S, MatrixX::rowi(Xa)) );			// treated as a column vector
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		for (; Xb!= Xend; ++Xb)
		{
			SymMatrix::value_type p = inner_prod(stemp, MatrixX::rowi(Xb));
			P(Xa.index1(),Xb.index1()) += p;
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
	typename MatrixX::const_iterator2 Xa = X.begin2();
	const typename MatrixX::const_iterator2 Xend = X.end2();
	typename MatrixX::const_iterator2 Xb;

	// P(a,b) = X'.row(a) * S * X.col(b) = X.col(a) * S * X.col(b)
	//        = X.col(b) * S * X.col(a)

	for (; Xa != Xend; ++Xa)				// Iterate Rows
	{
		stemp = ublas::prod(S, MatrixX::columni(Xa));			// treated as a column vector
		Xb = Xa;							// Start at the column Xa, only one triangle of symertric result required
		for (; Xb != Xend; ++Xb)
		{
			SymMatrix::value_type p = inner_prod(stemp, MatrixX::columni(Xb));
			P(Xa.index2(),Xb.index2()) += p;
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
	typename MatrixX::const_iterator2 Xa = X.begin2();
	const typename MatrixX::const_iterator2 Xend = X.end2();
	typename MatrixX::const_iterator2 Xb;

	// P(a,b) = X.col(a) * X.col(b)
	for (; Xa != Xend; ++Xa)				// Iterate vectors
	{		
		MatrixX::const_Column Xav = MatrixX::columni(Xa);
		Xb = Xa;							// Start at the row Xa only one triangle of symertric result required
		for (; Xb != Xend; ++Xb)
		{									// Tripple vector inner product
			SymMatrix::value_type p = 0;
			MatrixX::const_Column Xbv = MatrixX::columni(Xb);
			for (si = s.begin(); si != send; ++si) {	// TODO: use iterator on Xav and Xbv
				Vec::size_type i = si.index();
				p += Xav[i] * (*si) * Xbv[i];
			}
			P(Xa.index2(),Xb.index2()) += p;
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
	typename MatrixX::const_iterator2 Xa = X.begin2();
	const typename MatrixX::const_iterator2 Xend = X.end2();
	typename MatrixX::const_iterator2 Xb;

	// P(a,b) = X.row(a) * X.row(b)
	for (; Xa != Xend; ++Xa)				// Iterate vectors
	{		
		Xb = Xa;							// Start at the row Xa, only one triangle of symetric result required
		for (; Xb != Xend; ++Xb)
		{									// Simply multiple col Xa by col Xb
			SymMatrix::value_type p = inner_prod( MatrixX::columni(Xa), MatrixX::columni(Xb) );
			P(Xa.index2(),Xb.index2()) += p;
		}
	}
}


/*
 * prod_SPD - uBLAS Expression templates for X*X' and X*X'
 * Functions are only defined for the type for which the operation is efficient
 * ISSUE Although numerically symmetric, uBlas has no expression type to represent this property
 *  The result must be assigned to a symmetric container to exploit the symmetry
 */

template <class X>
struct prod_SPD_matrix_traits
{	// Provide ET result type XXT for prod(X,trans(X))
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_unary2_traits<X, ublas::scalar_identity<BOOST_UBLAS_TYPENAME X::value_type> >::result_type  XT_type;
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME X::value_type, X,
                                        BOOST_UBLAS_TYPENAME XT_type::value_type, XT_type>::result_type  XXT_type;

	// Provide ET result type XTX for prod(trans(X),X)
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME XT_type::value_type, XT_type,
                                        BOOST_UBLAS_TYPENAME X::value_type, X>::result_type  XTX_type;

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
 * ISSUE: Implemention only exploits sparseness in x
 * TODO: Define using ublas ET
 */
{
	Vec::value_type p = 0.;
	Vec::const_iterator xi = x.begin(), xi_end = x.end();
	for (; xi != xi_end; ++xi)
	{
		Vec::value_type x2 = *xi; x2 *= x2;
		p += x2 * s[xi.index()];
	}
	return p;
}


inline
prod_SPD_matrix_traits<ColMatrix>::XTX_type
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

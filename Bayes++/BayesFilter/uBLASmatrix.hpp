/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
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

#if (BOOST_VERSION >= 103100)
using ublas::noalias;
#else
namespace detail
{
/*
 * Improve syntax of effcient assignment where no aliases of LHS appear on the RHS
 *  noalias(lhs) = rhs_expression
 */
template <class C>
struct NoAliasAssign
{	
	NoAliasAssign (C& lval) : lval_(lval)
	{}
	template <class E>
	void operator= (const E& e)
	{	lval_.assign (e);
	}
	template <class E>
	void operator+= (const E& e)
	{	lval_.plus_assign (e);
	}
	template <class E>
	void operator-= (const E& e)
	{	lval_.minus_assign (e);
	}
	private:
	    C& lval_;
}; 
}//namespace detail

template <class E>
detail::NoAliasAssign<E> noalias(ublas::matrix_expression<E>& lvalue)
{
	return detail::NoAliasAssign<E> (lvalue() );
}
template <class E>
detail::NoAliasAssign<E> noalias(ublas::vector_expression<E>& lvalue)
{
	return detail::NoAliasAssign<E> (lvalue() );
}

#endif


namespace detail		// Lots of implementation detail
{

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
	explicit FMVec(const ublas::vector_expression<E>& e) : VecBase(e)
	{}	// vector_expression copy constructor
	template <class E>
	FMVec(const ublas::matrix_column<E>& e) : VecBase(e)
	{}	// vector_expression convsersion copy constructor, hides the implict copy required for matrix column access

	template <class E>
	FMVec& operator= (const ublas::vector_expression<E>& r)
	{	// Expression assignment; may be dependant on r
		VecBase::operator=(r);
		return *this;
	}
	FMVec& operator= (const FMVec& r)
	{	// Vector assignment; independant
		assign(r);
		return *this;
	}

	// Sub-range selection operators
	const ublas::vector_range<const VecBase> sub_range(size_t b, size_t e) const
	{
		return ublas::vector_range<const VecBase>(*this, ublas::range(b,e));
	}
	ublas::vector_range<VecBase> sub_range(size_t b, size_t e)
	{
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
	explicit FMMatrix(const ublas::matrix_expression<E>& e) : MatrixBase(e)
	{}	// matrix_expression constructor

	template <class E>
	FMMatrix& operator= (const ublas::matrix_expression<E>& r)
	{	// Expression assignment; may be dependant on r
		MatrixBase::operator=(r);
		return *this;
	}
	FMMatrix& operator= (const FMMatrix& r)
	{	// Matrix assignment; independant
		assign (r);
		return *this;
	}

	// Row,Column vector proxies
	typedef ublas::matrix_row<FMMatrix> Row;
	typedef const ublas::matrix_row<const FMMatrix> const_Row;
	typedef ublas::matrix_column<FMMatrix> Column;
	typedef const ublas::matrix_column<const FMMatrix> const_Column;

	// Vector proxies from iterators - static members dependant on MatrixBase type
	// ri() returns container associated with iterator. static_cast required as typeof(ri()) may not be typeof(MM)
	static Row rowi(const typename MatrixBase::iterator1& ri)
	{
		return Row(static_cast<FMMatrix&>(ri()), ri.index1());
	}
	static const_Row rowi(const typename MatrixBase::const_iterator1& ri)
	{
		return const_Row(static_cast<const FMMatrix&>(ri()), ri.index1());
	}
	static Column columni(const typename MatrixBase::iterator2& ci)
	{
		return Column(static_cast<FMMatrix&>(ci()), ci.index2());
	}
	static const_Column columni(const typename MatrixBase::const_iterator2& ci)
	{
		return const_Column(static_cast<const FMMatrix&>(ci()), ci.index2());
	}

	// Sub-range selection operators
	ublas::matrix_range<const MatrixBase>
	sub_matrix(size_t s1, size_t e1, size_t s2, size_t e2) const
	{
		return ublas::matrix_range<const MatrixBase> (*this, ublas::range(s1,e1), ublas::range(s2,e2));
	}
	ublas::matrix_range<MatrixBase>
	sub_matrix(size_t s1, size_t e1, size_t s2, size_t e2)
	{
		return ublas::matrix_range<MatrixBase> (*this, ublas::range(s1,e1), ublas::range(s2,e2));
	}

	// Requires boost_1.30.0 which has a generalised matrix_vector_slice
	ublas::matrix_vector_slice<const MatrixBase>
	sub_column(size_t s1, size_t e1, size_t s2) const 
	// Column vector s2 with rows [s1,e1)
	{
		return ublas::matrix_vector_slice<const MatrixBase> (*this, ublas::slice(s1,1,e1-s1), ublas::slice(s2,0,e1-s1));
	}
	ublas::matrix_vector_slice<MatrixBase>
	sub_column(size_t s1, size_t e1, size_t s2)
	// Column vector s2 with rows [s1,e1)
	{
		return ublas::matrix_vector_slice<MatrixBase> (*this, ublas::slice(s1,1,e1-s1), ublas::slice(s2,0,e1-s1));
	}
};


/*
 * Helper template to allow member construction before base class
 *  Boost version does not work as it passes by value
 */
template <typename MemberType>
class BaseFromMember
{
protected:
	MemberType member;
	explicit BaseFromMember() : member()
	{}

	template <typename T1>
	explicit BaseFromMember( const T1& x1 ) : member( x1 )
	{}

	template <typename T1, typename T2>
	explicit BaseFromMember( const T1& x1, const T2& x2 ) : member( x1, x2 )
	{}
};


/*
 * We require static type conversion between Symmetric matrices and equivilent row major matrices
 * Therefore we create symmetric matrix types, using a MatrixBase for storage
 * and wraps this in a symmetric_adaptor
 */
template <class MatrixBase>
class SymMatrixWrapper :
	private BaseFromMember<MatrixBase>,  // allow construction of rm before symmertic_adaptor
	public ublas::symmetric_adaptor<MatrixBase, ublas::upper>
{
	typedef BaseFromMember<MatrixBase> matrix_type;
	typedef ublas::symmetric_adaptor<MatrixBase, ublas::upper> adaptor_type;
public:
	SymMatrixWrapper () : matrix_type(), adaptor_type(member)
	{}
	SymMatrixWrapper (size_t nsize1, size_t nsize2) : matrix_type(nsize1,nsize2), adaptor_type(member)
	{}
	explicit SymMatrixWrapper (const SymMatrixWrapper& r) : matrix_type(reinterpret_cast<const MatrixBase&>(r)), adaptor_type(member)
	// Explict copy construction referencing the copy reinterpreted as a MatrixBase
	{}
	template <class E>
	explicit SymMatrixWrapper (const ublas::matrix_expression<E>& e) : matrix_type(e), adaptor_type(member)
	// Explict matrix_expression conversion constructor
	{}

	template <class E>
	SymMatrixWrapper& operator=(const ublas::matrix_expression<E>& r)
	{
		adaptor_type::operator=(r);
		return *this;
	}

	// Conversions straight to a FMMatrix, equivilent to a RowMatrix types
	const FMMatrix<MatrixBase>& asRowMatrix() const
	{
		return static_cast<const FMMatrix<MatrixBase>& >(member);
	}
	FMMatrix<MatrixBase>& asRowMatrix()
	{
		return static_cast<FMMatrix<MatrixBase>& >(member);
	}

	// Matrix storage members
	void clear()
	{	member.clear();
	}
	void resize(size_t nsize1, size_t nsize2)
	{
		member.resize(nsize1, nsize2);
	}
};

}//namespace detail




/*
 * Vector / Matrix types
 *  Finally the definitions !
 */
using detail::FMVec;		// Template class for template parameter matching
using detail::FMMatrix;

							// Default types
typedef FMVec<detail::BaseVector> Vec;
typedef FMMatrix<detail::BaseRowMatrix> RowMatrix;
typedef RowMatrix Matrix;
typedef FMMatrix<detail::BaseColMatrix> ColMatrix;
typedef FMMatrix<detail::SymMatrixWrapper<detail::BaseRowMatrix> > SymMatrix;
typedef FMMatrix<detail::BaseUpperTriMatrix> UTriMatrix;
typedef FMMatrix<detail::BaseLowerTriMatrix> LTriMatrix;
typedef FMMatrix<detail::BaseDiagMatrix> DiagMatrix;

							// Explicitly dense types
typedef FMVec<detail::BaseDenseVector> DenseVec;
typedef FMMatrix<detail::BaseDenseRowMatrix> DenseRowMatrix;
typedef DenseRowMatrix DenseMatrix;
typedef FMMatrix<detail::BaseDenseColMatrix> DenseColMatrix;
typedef FMMatrix<detail::SymMatrixWrapper<detail::BaseDenseRowMatrix> > DenseSymMatrix;
typedef FMMatrix<detail::BaseDenseUpperTriMatrix> DenseUTriMatrix;
typedef FMMatrix<detail::BaseDenseLowerTriMatrix> DenseLTriMatrix;
typedef FMMatrix<detail::BaseDenseDiagMatrix> DenseDiagMatrix;

							// Explicitly sparse types (any of the gappy types)
#ifdef BAYES_FILTER_GAPPY
typedef FMVec<detail::BaseSparseVector> SparseVec;
typedef FMMatrix<detail::BaseDenseRowMatrix> SparseRowMatrix;
typedef SparseRowMatrix SparseMatrix;
typedef FMMatrix<detail::BaseSparseColMatrix> SparseColMatrix;
typedef FMMatrix<detail::SymMatrixWrapper<detail::BaseSparseRowMatrix> > SparseSymMatrix;
#endif


/*
 * Matrix Adaptors, simply hide the uBLAS details
 */
template <class M>
const ublas::triangular_adaptor<const M, ublas::upper>
 UpperTri(const M& m)
/*
 * View Upper triangle of m
 * ISSUE VC7 cannot cope with UTriMatrix::functor1_type
 */
{
	return ublas::triangular_adaptor<const M, ublas::upper>(m);
}

template <class M>
const ublas::triangular_adaptor<const M, ublas::lower>
 LowerTri(const M& m)
/*
 * View Lower triangle of m
 */
{
	return ublas::triangular_adaptor<const M, ublas::lower>(m);
}



/*
 * Matrix Support Operations
 */
template <class Base>
ublas::matrix_vector_range<FMMatrix<Base> >
 diag(FMMatrix<Base>& M, size_t n)
{	// Return a vector proxy to the first n diagonal elements of M
	return ublas::matrix_vector_range<FMMatrix<Base> >(M, ublas::range(0,n), ublas::range(0,n));
}

template <class Base>
const ublas::matrix_vector_range<const FMMatrix<Base> >
 diag(const FMMatrix<Base>& M, size_t n)
{	// Return a vector proxy to the first n diagonal elements of M
	return ublas::matrix_vector_range<const FMMatrix<Base> >(M, ublas::range(0,n), ublas::range(0,n));
}

template <class Base>
ublas::matrix_vector_range<FMMatrix<Base> >
 diag(FMMatrix<Base>& M)
{	// Return a vector proxy to the diagonal elements of M
	const size_t common_size = std::min(M.size1(),M.size2());
	return ublas::matrix_vector_range<FMMatrix<Base> >(M, ublas::range(0,common_size), ublas::range(0,common_size));
}

template <class Base>
const ublas::matrix_vector_range<const FMMatrix<Base> >
 diag(const FMMatrix<Base>& M)
{	// Return a vector proxy to the diagonal elements of M
	const size_t common_size = std::min(M.size1(),M.size2());
	return ublas::matrix_vector_range<const FMMatrix<Base> >(M, ublas::range(0,common_size), ublas::range(0,common_size));
}

template <class Base>
void identity(FMMatrix<Base>& I)
{	// Set I to generalised Identity matrix. Clear I and set diag(I) to one
	I.clear();
							// Set common diagonal elements
	size_t common_size = std::min(I.size1(),I.size2());
	typedef typename Base::value_type Base_value_type;
	diag(I) = ublas::scalar_vector<Base_value_type>(common_size, Base_value_type(1));
}



/*
 * Symmetric Positive (Semi) Definate multiplication: X*S*X' and X'*S*X
 *  The name is slightly misleading. The result is actually only PD if S is PD
 *  Algorithms are intended to exploit the symmerty of the result
 *  and also where possible the row by row multiplication inherent in X*X'
 */

namespace detail {		// mult_SPD now an implementation detail
	
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
		typename MatrixX::const_Row Xav = MatrixX::rowi(Xa);
		Xb = Xa;							// Start at the row Xa only one triangle of symetric result required
		for (; Xb != Xend; ++Xb)
		{
			SymMatrix::value_type p = 0;	// Tripple vector inner product
			typename MatrixX::const_Row Xbv = MatrixX::rowi(Xb);
			for (si = s.begin(); si != send; ++si) {
				Vec::size_type i = si.index();
				p += Xav[i] * (*si) * Xbv[i];
			}
			P(Xa.index1(),Xb.index1()) += p;
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
		typename MatrixX::const_Column Xav = MatrixX::columni(Xa);
		Xb = Xa;							// Start at the row Xa only one triangle of symertric result required
		for (; Xb != Xend; ++Xb)
		{									// Tripple vector inner product
			SymMatrix::value_type p = 0;
			typename MatrixX::const_Column Xbv = MatrixX::columni(Xb);
			for (si = s.begin(); si != send; ++si) {
				Vec::size_type i = si.index();
				p += Xav[i] * (*si) * Xbv[i];
			}
			P(Xa.index2(),Xb.index2()) += p;
		}
	}
}

}//namespace detail


/*
 * prod_SPD - uBLAS Expression templates for X*X' and X*X'
 * Functions are only defined for the type for which the operation is efficient
 * ISSUE Although numerically symmetric, uBlas has no expression type to represent this property
 *  The result must be assigned to a symmetric container to exploit the symmetry
 */

template <class E1, class E2>
struct prod_expression_result
{	// Provide ET result E1E2T_type of prod(matrix_expression<E1>,trans(matrix_expression<E2>)
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_unary2_traits<E2, ublas::scalar_identity<BOOST_UBLAS_TYPENAME E2::value_type> >::result_type  E2T_type;
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                        BOOST_UBLAS_TYPENAME E2T_type::value_type, E2T_type>::result_type  E1E2T_type;

	// Provide ET result E1TE2_type of prod(trans(matrix_expression<E1>),matrix_expression<E2>)
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_unary2_traits<E1, ublas::scalar_identity<BOOST_UBLAS_TYPENAME E1::value_type> >::result_type  E1T_type;
	typedef BOOST_UBLAS_TYPENAME ublas::matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME E1T_type::value_type, E1T_type,
                                        BOOST_UBLAS_TYPENAME E2::value_type, E2>::result_type  E1TE2_type;
};

 
template<class E> inline
typename prod_expression_result<E,E>::E1E2T_type
 prod_SPD (const ublas::matrix_expression<E>& X)
/*
 * Symmetric Positive (Semi) Definate product: X*X'
 */
{
	typedef typename E::value_type t1;
	return prod( X, trans(X) );
}

template<class E1, class E2> inline
typename prod_expression_result<E1,E2>::E1E2T_type
 prod_SPD (const ublas::matrix_expression<E1>& X, const SymMatrix& S, ublas::matrix_expression<E2>& XStemp)
/*
 * Symmetric Positive (Semi) Definate product: X*(X*S)', XStemp = X*S
 */
{
	return prod( X, trans(XStemp().assign(prod(X,S))) );
}

inline
SymMatrix prod_SPD (const RowMatrix& X, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate product: X*diag_matrix(s)*X'
 * TODO: Define using ublas ET
 */
{
	SymMatrix P(X.size1(),X.size1());
	P.clear();
	detail::mult_SPD(X,s, P);
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


template<class E> inline
typename prod_expression_result<E,E>::E1TE2_type
 prod_SPDT (const ublas::matrix_expression<E>& X)
/*
 * Symmetric Positive (Semi) Definate product: X'*X
 */
{
	return prod( trans(X), X );
}

inline
SymMatrix prod_SPDT (const ColMatrix& X, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate product: X'*diag_matrix(s)*X
 * TODO: Define using ublas ET
 */
{
	SymMatrix P(X.size2(),X.size2());
	P.clear();
	detail::mult_SPDT(X,s, P);
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

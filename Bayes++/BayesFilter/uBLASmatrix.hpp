/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
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
 *  The computation of symetric products P*P' is can be efficiently computed with row operations. 
 *  These products are used extensively in covariance calculation. Therefore we use
 *  row_major as the default matrix representation.
 */

/* Filter Matrix Namespace */
namespace Bayesian_filter_matrix
{
						// Allow use a few functions in own namespace (particular useful for compilers with Konig lookup)
using ublas::row;
using ublas::column;
using ublas::trans;
using ublas::prod;		// These do not apply to the templated prod<temp> funtions
using ublas::inner_prod;
using ublas::outer_prod;
using ublas::noalias;

enum EmptyTag {Empty};	// Tag type used for empty matrix constructor


namespace detail		// Lots of implementation detail
{

/*
 * Filter Vec type
 */
template <class VecBase>
class FMVec : public VecBase
{
public:
	typedef typename VecBase::value_type value_type;

	// No Default Constructor. Empty creation is very error prone
	explicit FMVec(EmptyTag) : VecBase()
	{}	// Empty constructor
	explicit FMVec(typename VecBase::size_type size) : VecBase(size)
	{	// Sized constructor
		assign (ublas::scalar_vector<value_type>(size, std::numeric_limits<value_type>::signaling_NaN()));
	}
	FMVec(const FMVec& c) : VecBase(static_cast<const VecBase&>(c))
	{}	// Copy constructor
	template <class E>
	explicit FMVec(const ublas::vector_expression<E>& e) : VecBase(e)
	{}	// vector_expression copy constructor
	template <class E>
	FMVec(const ublas::matrix_column<E>& e) : VecBase(e)
	{}	// conversion copy constructor, hides the implict copy required for matrix column access

	template <class E>
	FMVec& operator= (const ublas::vector_expression<E>& r)
	{	// Expression assignment; may be dependant on r
		VecBase::operator=(r);
		return *this;
	}
	FMVec& operator= (const FMVec& r)
	{	// Vector assignment; independant
		VecBase::operator=(static_cast<const VecBase&>(r));
		return *this;
	}

	// Sub-range selection operators
	const ublas::vector_range<const VecBase> sub_range(typename VecBase::size_type b, typename VecBase::size_type e) const
	{
		return ublas::vector_range<const VecBase>(*this, ublas::range(b,e));
	}
	ublas::vector_range<VecBase> sub_range(typename VecBase::size_type b, typename VecBase::size_type e)
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
	FMMatrix(typename MatrixBase::size_type size1, typename MatrixBase::size_type size2) : MatrixBase(size1,size2)
	{	// Sized constructor
		assign (ublas::scalar_matrix<value_type>(size1, size2, std::numeric_limits<value_type>::signaling_NaN()));
	}
	FMMatrix(const FMMatrix& c) : MatrixBase(static_cast<const MatrixBase&>(c))
	{}	// Copy constructor
	template <class E>
	explicit FMMatrix(const ublas::matrix_expression<E>& e) : MatrixBase(e)
	{}	// matrix_expression copy constructor

	template <class E>
	FMMatrix& operator= (const ublas::matrix_expression<E>& r)
	{	// Expression assignment; may be dependant on r
		MatrixBase::operator=(r);
		return *this;
	}
	FMMatrix& operator= (const FMMatrix& r)
	{	// Matrix assignment; independant
		MatrixBase::operator=(static_cast<const MatrixBase&>(r));
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
	sub_matrix(typename MatrixBase::size_type s1, typename MatrixBase::size_type e1, typename MatrixBase::size_type s2, typename MatrixBase::size_type e2) const
	{
		return ublas::matrix_range<const MatrixBase> (*this, ublas::range(s1,e1), ublas::range(s2,e2));
	}
	ublas::matrix_range<MatrixBase>
	sub_matrix(typename MatrixBase::size_type s1, typename MatrixBase::size_type e1, typename MatrixBase::size_type s2, typename MatrixBase::size_type e2)
	{
		return ublas::matrix_range<MatrixBase> (*this, ublas::range(s1,e1), ublas::range(s2,e2));
	}

	// Requires boost_1.30.0 which has a generalised matrix_vector_slice
	ublas::matrix_vector_slice<const MatrixBase>
	sub_column(typename MatrixBase::size_type s1, typename MatrixBase::size_type e1, typename MatrixBase::size_type s2) const 
	// Column vector s2 with rows [s1,e1)
	{
		return ublas::matrix_vector_slice<const MatrixBase> (*this, ublas::slice(s1,1,e1-s1), ublas::slice(s2,0,e1-s1));
	}
	ublas::matrix_vector_slice<MatrixBase>
	sub_column(typename MatrixBase::size_type s1, typename MatrixBase::size_type e1, typename MatrixBase::size_type s2)
	// Column vector s2 with rows [s1,e1)
	{
		return ublas::matrix_vector_slice<MatrixBase> (*this, ublas::slice(s1,1,e1-s1), ublas::slice(s2,0,e1-s1));
	}
};


/*
 * Helper template to allow member construction before base class
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
	private BaseFromMember<MatrixBase>,  // allow construction of MatrixBase member before symmetric_adaptor
	public ublas::symmetric_adaptor<MatrixBase, ublas::upper>
{
	typedef BaseFromMember<MatrixBase> matrix_type;
	typedef ublas::symmetric_adaptor<MatrixBase, ublas::upper> symadaptor_type;
public:
	SymMatrixWrapper () : matrix_type(), symadaptor_type(matrix_type::member)
	{}
	SymMatrixWrapper (typename MatrixBase::size_type size1, typename MatrixBase::size_type size2) : matrix_type(size1,size2), symadaptor_type(matrix_type::member)
	{}	// Normal sized constructor
	explicit SymMatrixWrapper (const SymMatrixWrapper& r) : matrix_type(reinterpret_cast<const MatrixBase&>(r)), symadaptor_type(matrix_type::member)
	{}	// Explict copy construction referencing the copy reinterpreted as a MatrixBase
	template <class E>
	explicit SymMatrixWrapper (const ublas::matrix_expression<E>& e) : matrix_type(e), symadaptor_type(matrix_type::member)
	{}	// Explict matrix_expression conversion constructor

	template <class E>
	SymMatrixWrapper& operator=(const ublas::matrix_expression<E>& r)
	{
		symadaptor_type::operator=(r);
		return *this;
	}

	// Conversions straight to a FMMatrix, equivilent to a RowMatrix types
	const FMMatrix<MatrixBase>& asRowMatrix() const
	{
		return static_cast<const FMMatrix<MatrixBase>& >(matrix_type::member);
	}
	FMMatrix<MatrixBase>& asRowMatrix()
	{
		return static_cast<FMMatrix<MatrixBase>& >(matrix_type::member);
	}

	// Matrix storage members
	void clear()
	{
		matrix_type::member.clear();
	}
	void resize(typename MatrixBase::size_type nsize1, typename MatrixBase::size_type nsize2)
	{
		matrix_type::member.resize(nsize1, nsize2);
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
 diag(FMMatrix<Base>& M, typename Base::size_type n)
{	// Return a vector proxy to the first n diagonal elements of M
	return ublas::matrix_vector_range<FMMatrix<Base> >(M, ublas::range(0,n), ublas::range(0,n));
}

template <class Base>
const ublas::matrix_vector_range<const FMMatrix<Base> >
 diag(const FMMatrix<Base>& M, typename Base::size_type n)
{	// Return a vector proxy to the first n diagonal elements of M
	return ublas::matrix_vector_range<const FMMatrix<Base> >(M, ublas::range(0,n), ublas::range(0,n));
}

template <class Base>
ublas::matrix_vector_range<FMMatrix<Base> >
 diag(FMMatrix<Base>& M)
{	// Return a vector proxy to the diagonal elements of M
	const typename Base::size_type common_size = std::min(M.size1(),M.size2());
	return ublas::matrix_vector_range<FMMatrix<Base> >(M, ublas::range(0,common_size), ublas::range(0,common_size));
}

template <class Base>
const ublas::matrix_vector_range<const FMMatrix<Base> >
 diag(const FMMatrix<Base>& M)
{	// Return a vector proxy to the diagonal elements of M
	const typename Base::size_type common_size = std::min(M.size1(),M.size2());
	return ublas::matrix_vector_range<const FMMatrix<Base> >(M, ublas::range(0,common_size), ublas::range(0,common_size));
}

template <class Base>
void identity(FMMatrix<Base>& I)
{	// Set I to generalised Identity matrix. Clear I and set diag(I) to one
	I.clear();
							// Set common diagonal elements
	typename Base::size_type common_size = std::min(I.size1(),I.size2());
	typedef typename Base::value_type Base_value_type;
	diag(I).assign (ublas::scalar_vector<Base_value_type>(common_size, Base_value_type(1)));
}


/*
 * prod_SPD - uBLAS Expression templates for X*X' , X*S*X´ , X'*X and 'X*S*X
 *  The name is slightly misleading. The result is actually only PD if S is PD
 * Algorithms are intended to exploit the symmerty of the result and S, therefore
 * are often defined with a more efficient compuation order.
 * Functions are only defined for the type for which the operation is efficient
 * ISSUE Although numerically symmetric, uBlas has no expression type to represent this property
 *  The result must be assigned to a symmetric container to exploit the symmetry
 */

template <class E1, class E2>
struct prod_expression_result
{	// Provide ET result E1E2T_type of prod(matrix_expression<E1>,trans(matrix_expression<E2>)
	typedef typename ublas::matrix_unary2_traits<E2, ublas::scalar_identity<typename E2::value_type> >::result_type  E2T_type;
	typedef typename ublas::matrix_matrix_binary_traits<typename E1::value_type, E1,
                                        typename E2T_type::value_type, E2T_type>::result_type  E1E2T_type;

	// Provide ET result E1TE2_type of prod(trans(matrix_expression<E1>),matrix_expression<E2>)
	typedef typename ublas::matrix_unary2_traits<E1, ublas::scalar_identity<typename E1::value_type> >::result_type  E1T_type;
	typedef typename ublas::matrix_matrix_binary_traits<typename E1T_type::value_type, E1T_type,
                                        typename E2::value_type, E2>::result_type  E1TE2_type;
};

 
template<class E> inline
typename prod_expression_result<E,E>::E1E2T_type
 prod_SPD (const ublas::matrix_expression<E>& X)
/*
 * Symmetric Positive (Semi) Definate product: X*X'
 */
{
	// ISSUE: uBLAS post Boost 1_31_0 introduces a trans(const matrix_expression<E>& e) which propogates the const expression type
	// Bypass this to avoid having to detect the Boost version
	return prod( X, trans(const_cast<ublas::matrix_expression<E>&>(X) ));
}

template<class EX, class ES, class ET> inline
typename prod_expression_result<EX,ET>::E1E2T_type
 prod_SPD (const ublas::matrix_expression<EX>& X, const ublas::matrix_expression<ES>& S, ublas::matrix_expression<ET>& XStemp)
/*
 * Symmetric Positive (Semi) Definate product: X*(X*S)', XStemp = X*S
 * Assumes S symmetric
 */
{
	return prod( X, trans(prod(X,S,XStemp())) );
}

template<class EX, class ES, class ET> inline
ET prod_SPD (const ublas::matrix_expression<EX>& X, const ublas::vector_expression<ES>& s, ublas::matrix_expression<ET>& Ptemp)
/*
 * Symmetric Positive (Semi) Definate product: X*diag_matrix(s)*X', Ptemp = return value
 * Precond: Ptemp must be size conformant with the product
 * TODO requires a prod_diag(X,s)
 */
{
	const EX& XX = X();
	Vec::const_iterator si, send = s().end();
	typename EX::const_iterator1 Xa = XX.begin1();
	const typename EX::const_iterator1 Xend = XX.end1();
	typename EX::const_iterator1 Xb;

	// P(a,b) = sum(X.row(a) * s * X.row(b))
	for (; Xa != Xend; ++Xa)		// Iterate rows
	{
		typedef const ublas::matrix_row<const EX> EX_row;
		EX_row Xav (ublas::row_const(XX, Xa.index1()));
		Xb = Xa;							// Start at the row Xa only one triangle of symetric result required
		for (; Xb != Xend; ++Xb)
		{
			typename EX::value_type p = 0;	// Triple vector inner product
			EX_row Xbv (ublas::row_const(XX, Xb.index1()));
			for (si = s().begin(); si != send; ++si) {
				Vec::size_type i = si.index();
				p += Xav[i] * (*si) * Xbv[i];
			}
			Ptemp()(Xa.index1(),Xb.index1()) = p;
			Ptemp()(Xb.index1(),Xa.index1()) = p;
		}
	}
	return Ptemp();
}


template<class E> inline
typename prod_expression_result<E,E>::E1TE2_type
 prod_SPDT (const ublas::matrix_expression<E>& X)
/*
 * Symmetric Positive (Semi) Definate product: X'*X
 */
{
	// ISSUE: See prod_SPD
	return prod( trans(const_cast<ublas::matrix_expression<E>&>(X) ), X);
}

template<class EX, class ES, class ET> inline
typename prod_expression_result<ET,EX>::E1TE2_type
 prod_SPDT (const ublas::matrix_expression<EX>& X, const ublas::matrix_expression<ES>& S, ublas::matrix_expression<ET>& SXtemp)
/*
 * Symmetric Positive (Semi) Definate product: (S*X)'*X, SXtemp = S*X
 * Assumes S symmetric
 */
{
	return prod( trans(prod(S,X,SXtemp())), X);
}

template<class EX, class ES, class ET> inline
ET prod_SPDT (const ublas::matrix_expression<EX>& X, const ublas::vector_expression<ES>& s, ublas::matrix_expression<ET>& Ptemp)
/*
 * Symmetric Positive (Semi) Definate product: X'*diag_matrix(s)*X, Ptemp = return value
 * Precond: Ptemp must be size conformant with the product
 * TODO requires a prod_diag(X,s)
 */
{
	const EX& XX = X();
	Vec::const_iterator si, send = s().end();
	typename EX::const_iterator2 Xa = X.begin2();
	const typename EX::const_iterator2 Xend = X.end2();
	typename EX::const_iterator2 Xb;

	// P(a,b) = sum(X.col(a) * s * X.col(b))
	for (; Xa != Xend; ++Xa)		// Iterate columns
	{
		typedef const ublas::matrix_column<const EX> EX_column;
		EX_column Xav = ublas::column_const(XX, Xa.index2());
		Xb = Xa;							// Start at the column Xa only one triangle of symetric result required
		for (; Xb != Xend; ++Xb)
		{
			typename EX::value_type p = 0;	// Triple vector inner product
			EX_column Xbv = ublas::column_const(XX, Xb.index2());
			for (si = s().begin(); si != send; ++si) {
				Vec::size_type i = si.index();
				p += Xav[i] * (*si) * Xbv[i];
			}
			Ptemp()(Xa.index2(),Xb.index2()) = p;
			Ptemp()(Xb.index2(),Xa.index2()) = p;
		}
	}
	return Ptemp();
}

inline Vec::value_type prod_SPDT (const Vec& x, const Vec& s)
/*
 * Symmetric Positive (Semi) Definate product: x'*diag_matrix(s)*x
 */
{
	return inner_prod(x, ublas::element_prod(s,x));
}

inline Vec::value_type prod_SPDT (const Vec& x, const SymMatrix& S)
/*
 * Symmetric Positive (Semi) Definate multiply:	p = x'*S*x
 */
{
	return inner_prod (x, prod(S,x) );
}

}//namespace

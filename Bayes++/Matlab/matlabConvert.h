/*
 * Bayesian_filter_matrix <-> Matlab conversions
 */

namespace Matlab_convert
{

template <class FilterMatrix>
void Array(mxArray* ma, const FilterMatrix& M)
/*
 * Convert Filter Matrix into a mxArray
 */
{
	FM::Subscript rows = M.size1(), cols = M.size2();

	double *mai = mxGetPr(ma);
	for (FM::Subscript c = 0; c < cols; ++c)
	{
		for (FM::Subscript r = 0; r < rows; ++r)
		{
			*mai++ = M(r,c);
		}
	}
}

template <class FilterMatrix>
mxArray* Array(const FilterMatrix& M)
/*
 * Convert Filter Matrix into a mxArray
 */
{
	FM::Subscript rows = M.size1(), cols = M.size2();
	mxArray* ma = mxCreateDoubleMatrix(rows,cols, mxREAL);

	Array (ma, M);
	return ma;
}

template <class FilterMatrix>
mxArray* ArrayTranspose(const FilterMatrix& M)
/*
 * Convert Filter Matrix's transpose into a mxArray
 */
{
	FM::Subscript rows = M.size1(), cols = M.size2();
	mxArray* ma = mxCreateDoubleMatrix(cols,rows, mxREAL);

	double *mai = mxGetPr(ma);
	for (FM::Subscript r = 0; r < rows; ++r)
	{
		for (FM::Subscript c = 0; c < cols; ++c)
		{
			*mai++ = M(r,c);
		}
	}
	return ma;
}

void Array(mxArray* ma, const Bayesian_filter_matrix::Vec& V)
/*
 * Convert Vec into a mxArray column vector
 */
{
	FM::Subscript size = V.size();
	double *mai = mxGetPr(ma);
	for (FM::Subscript r = 0; r < size; ++r)
	{
		*mai++ = V[r];
	}
}

mxArray* Array(const Bayesian_filter_matrix::Vec& V)
/*
 * Convert Vec into a mxArray column vector
 */
{
	FM::Subscript size = V.size();
	mxArray* ma = mxCreateDoubleMatrix(size,1, mxREAL);

	Array (ma, V);
	return ma;
}

mxArray* Array(const double d)
/*
 * Convert double into a mxArray scalar
 */
{
	mxArray* ma = mxCreateDoubleMatrix(1,1, mxREAL);

	double *mai = mxGetPr(ma);
	*mai = d;
	return ma;
}

mxArray* Array(const unsigned u)
/*
 * Convert unsigned into a mxArray scalar
 */
{
	mxArray* ma = mxCreateDoubleMatrix(1,1, mxREAL);

	double *mai = mxGetPr(ma);
	*mai = double(u);
	return ma;
}

Bayesian_filter_matrix::Vec
Vector(const mxArray* ma)
/*
 * Convert mxArray a Filter Vector
 */
{
				// Input must be a row/column vector of double
	int rows = mxGetM(ma);
	int cols = mxGetN(ma);
    if (rows != 1 && cols != 1)
		mexErrMsgTxt("Expecting a row/column vector");
	if (!mxIsDouble(ma) || mxIsComplex(ma))
		mexErrMsgTxt("Vector must be of non-complex scalar double");
    
	int size = rows > cols ? rows : cols;
	Bayesian_filter_matrix::Vec V(size);

	double *mai = mxGetPr(ma);
	for (int i = 0; i < size; ++i)
	{
		V[i] = *mai++;
	}
	return V;
};

Bayesian_filter_matrix::Matrix
Matrix(const mxArray* ma)
/*
 * Convert mxArray to a Filter Matrix
 */
{
				// Input must be an array of double
	if (!mxIsDouble(ma) || mxIsComplex(ma))
		mexErrMsgTxt("Array must be of non-complex scalar double");

	int rows = mxGetM(ma);
	int cols = mxGetN(ma);
	Bayesian_filter_matrix::Matrix M(rows,cols);

	double *mai = mxGetPr(ma);
	for (int c = 0; c < cols; ++c)
	{
		for (int r = 0; r < rows; ++r)
		{
			M(r,c) = *mai++;
		}
	}
	return M;
};

Bayesian_filter_matrix::SymMatrix
SymMatrix(const mxArray* ma)
/*
 * Convert mxArray to a Filter SymMatrix
 */
{
				// Input must be an array of double
	if (!mxIsDouble(ma) || mxIsComplex(ma))
		mexErrMsgTxt("Array must be of non-complex scalar double");

	int rows = mxGetM(ma);
	int cols = mxGetN(ma);
	if (rows != cols)
	{
		mexErrMsgTxt("Array must be square");
	}

	Bayesian_filter_matrix::SymMatrix M(rows,cols);

	double *mai = mxGetPr(ma);
	for (int c = 0; c < cols; ++c)
	{
		int r;
		for (r = 0; r < c; ++r)
			*mai++;		// Ignore lower triangle
		for (;r < rows; ++r)
		{
			M(r,c) = *mai++;
		}
	}
	return M;
};

char* cstring (const mxArray* ma)
/*
 * Convert mxArray to null terminated C string
  * Return NULL for non string type
 */
{
    /* must be a string. */
	if (mxIsChar(ma) != 1)
		return NULL;
    /* Input must be a row vector. */
    if (mxGetM(ma) != 1)
		return NULL;      
    
    /* Get the length of the input string. */
    int numchar = (mxGetM(ma) * mxGetN(ma)) + 1;

	/* Allocate memory for input and output strings: no error return in MEX */
    char* buf=(char*)mxCalloc(numchar, sizeof(mxChar));

    /* Copy the string data into a C string */
	int status = mxGetString(ma, buf, numchar*sizeof(mxChar));
	if (status)
		mexErrMsgTxt("Internal string conversion error");
	return buf;
}

double cdouble (const mxArray* ma)
/*
 * Convert mxArray to double
 * Error if mxArray is not of expected type
 */
{
	/* The input must be a noncomplex scalar double.*/
	if( !mxIsDouble(ma) || mxIsComplex(ma) || !(mxGetM(ma)==1 && mxGetN(ma)==1) )
	{
		mexErrMsgTxt("Expecting a noncomplex scalar double.");
	}
	return *mxGetPr(ma);
}  

inline int cint (const mxArray* ma)
/*
 * Convert mxArray to int
 * Error if mxArray is not of expected type. This is actuall a noncomplex scalar double
 * as there is no way for matlab to create an int array other then using MEX
 */
{
	return int(cdouble(ma));
}

inline int cunsigned (const mxArray* ma)
/*
 * Convert mxArray to unsigned
 * Error if mxArray is not of expected type. This is actuall a noncomplex scalar double
 * as there is no way for matlab to create an int array other then using MEX
 */
{
	return unsigned(cdouble(ma));
}

}//namespace Matlab_convert

/*
 * Bayesian_filter_matrix <-> MuPAD conversions
 *  Many MuPAD function use long, we assume long === int
 */

namespace MuPAD_convert
{

template <class FilterMatrix>
MTcell List(const FilterMatrix& M)
/*
 * Convert Filter Matrix into a MuPAD list
 */
{
	MTcell mpList = MFnewList(M.size1());
    // creates a MuPAD list of lists
	for (unsigned r=0; r < M.size1(); ++r)
	{
	    // creates a MuPAD list
		MTcell colList = MFnewList(M.size2());
		for (unsigned c=0; c < M.size2(); ++c) {
			MFsetList( &colList, c, MFdouble(M(r,c)));
		}
		MFsig(colList);
		MFsetList (&mpList, r, colList);
	}
	MFsig (mpList);                         // calculates the signature of the list
	return mpList;
}

template <>
MTcell List(const Bayesian_filter_matrix::Vec& V)
/*
 * Convert Vector into a MuPAD list
 */
{
	MTcell mpList = MFnewList(V.size());
    // creates a MuPAD list
	for (unsigned r=0; r < V.size(); ++r)
	{
		MFsetList (&mpList, r, MFdouble(V[r]) );
	}
	MFsig (mpList);                         // calculates the signature of the list
	return mpList;
}

template <class FilterMatrix>
MTcell ListTranspose(const FilterMatrix& M)
/*
 * Convert Filter Matrix's transpose into a MuPAD list
 */
{
	MTcell mpList = MFnewList(M.size2());
    // creates a MuPAD list of lists
	for (unsigned r=0; r < M.size2(); ++r)
	{
	    // creates a MuPAD list
		MTcell colList = MFnewList(M.size1());
		for (unsigned c=0; c < M.size1(); ++c) {
			MFsetList( &colList, c, MFdouble(M(c,r)));
		}
		MFsig(colList);
		MFsetList (&mpList, r, colList);
	}
	MFsig (mpList);                         // calculates the signature of the list
	return mpList;
}


template <class FilterMatrix>
MTcell Array(const FilterMatrix& M)
/*
 * Convert Filter Matrix into a MuPAD Array
 */
{
	MTcell mpList = List(M);
	// Convert to array
	MTcell mpArray = MFlist2array(mpList, M.size1(),M.size2());
	MFfree(mpList);
	return mpArray;
}

template <>
MTcell Array(const Bayesian_filter_matrix::Vec& V)
/*
 * Convert Vector into a MuPAD Array
 */
{
	MTcell mpList = List(V);
	// Convert to array
	MTcell mpArray = MFlist2array(mpList, V.size());
	MFfree(mpList);
	return mpArray;
}

template <class FilterMatrix>
MTcell ArrayTranspose(const FilterMatrix& M)
/*
 * Convert Filter Matrix's transpose into a MuPAD Array
 */
{
	MTcell mpList = ListTranspose(M);
	// Convert to array
	MTcell mpArray = MFlist2array(mpList, M.size2(),M.size1());
	MFfree(mpList);
	return mpArray;
}


template <class Dummy>
Bayesian_filter_matrix::Vec
Vector(Dummy);

template <>
Bayesian_filter_matrix::Vec
Vector(MTcell aorl)
/*
 * Convert MuPAD list OR array to a Filter Vector
 */
{
	// Check type
	if (MFisList(aorl))
	{
		// list is simple
		int size = MFnops(aorl);
		Bayesian_filter_matrix::Vec V(size);

		for (int i = 0; i < size; ++i)
		{
			V[i] = MFdouble (MFgetList(&aorl,i));
		}
		return V;
	}
	else
	{
		// not a list, accept 1 dim arrays
		if (!MFisArray(aorl))
			MFerror( "array or list expected" );
		int dim = MFdimArray(aorl);
		if (dim == 1 || dim == 2)
		{
			// 1d row Array or 2d col Array
			long lb, ub;
			MFrangeArray(aorl, 1, &lb, &ub);
			int size = (ub-lb+1);
			if (dim == 2)
			{
				MFrangeArray(aorl, 2, &lb, &ub);
				if ((ub-lb+1) != 1)
					MFerror( "2d array with 1 column required for vector" );
			}

			MTcell list = MFarray2list(aorl);

			// extract list into Matrix
			Bayesian_filter_matrix::Vec V(size);
			for (int j = 0; j < size; ++j)
			{
				V[j] = MFdouble (MFgetList(&list,j));
			}

			MFfree (list);							// free list
			return V;
		}
		else {
			MFerror( "1d or 2d array required for vector" );
			// Unreachable
			Bayesian_filter_matrix::Vec V(Bayesian_filter_matrix::Empty);
			return V;
		}
	}
};

template <class Dummy>
Bayesian_filter_matrix::Matrix
Matrix_from_list(Dummy);

template <>
Bayesian_filter_matrix::Matrix
Matrix_from_list(MTcell list)
/*
 * Convert MuPAD list to a Filter Matrix
 */
{
	bool conversionFailed = false;
	
	// find size and check all sub-lists of same size?
	int rows = MFnops(list);
	int cols;
	if (MFisList(list))
	{
		if (rows == 0) {
			cols = 0;
		}
		else
		{
			cols = MFnops( MFgetList(&list,0) );
			for (int i=0; i < rows && !conversionFailed; ++i)
			{
				if (!MFisList(MFgetList(&list,i)) || 
					!MFnops  (MFgetList(&list,i)) == cols)
				{
					conversionFailed = true;
				}
			}
		}
	}
	else
		conversionFailed = true;

	if (conversionFailed)
		MFerror( "Conversion of list of list to Matrix failed" );

	Bayesian_filter_matrix::Matrix M(rows, cols);
	for (int i=0; i < rows; ++i) 
	{
		MTcell sublist = MFgetList(&list,i);
		for (int j = 0; j < cols; ++j)
		{
			M(i,j) = MFdouble (MFgetList(&sublist,j));
		}
	}
	
	
	return M;
}

template <class Dummy>
Bayesian_filter_matrix::Matrix
Matrix(Dummy);

template <>
Bayesian_filter_matrix::Matrix
Matrix(MTcell aorl)
/*
 * Convert MuPAD list OR array to a Filter Matrix
 */
{
	// Check type
	if (MFisList(aorl))
	{
		// list is simple
		Bayesian_filter_matrix::Matrix& M = Matrix_from_list(aorl);
		return M;
	}
	else
	{
		// not a list, accept 2 dim arrays
		if (!MFisArray(aorl))
			MFerror( "array or list expected" );
		int dim = MFdimArray(aorl);
		if (dim != 2)
			MFerror( "2d array expected" );

		long lb, ub;
		MFrangeArray(aorl, 1, &lb, &ub);
		int cols = (ub-lb+1);
		MFrangeArray(aorl, 2, &lb, &ub);
		int rows = (ub-lb+1);

		// convert into list
		MTcell list = MFarray2list(aorl);

		// extract list into Matrix
		Bayesian_filter_matrix::Matrix M(rows, cols);
		int listit = 0;
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				M(i,j) = MFdouble (MFgetList(&list,listit++));
			}

		}
		MFfree (list);							// free list

		return M;
	}
};


}//namespace MuPAD_convert

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 *  Kalman filter representing representation of SLAM
 */

		// Bayes++ Bayesian filtering schemes
#include "BayesFilter/bayesFlt.hpp"
		// Bayes++ SLAM
#include "SLAM.hpp"
#include "kalmanSLAM.hpp"

namespace SLAM_filter
{

template <class Base>
inline void zero(FM::ublas::matrix_range<Base> A)
// Zero a matrix_range
{	// Note A cannot be a reference
	typedef typename Base::value_type Base_value_type;
	FM::noalias(A) = FM::ublas::scalar_matrix<Base_value_type>(A.size1(),A.size2(), Base_value_type());
}

Kalman_SLAM::Kalman_SLAM( Kalman_filter_generator& filter_generator ) :
	SLAM(),
	fgenerator(filter_generator),
	full(NULL)
{
	nL = 0;
	nM = 0;
}

Kalman_SLAM::~Kalman_SLAM()
{
	fgenerator.dispose (full);
}

void Kalman_SLAM::init_kalman (const FM::Vec& x, const FM::SymMatrix& X)
{
	// TODO maintain map states on reinit
	nL = x.size();
	nM = 0;
	if (full) fgenerator.dispose (full);
		// generate full filter
	full = fgenerator.generate(nL);
		// initialise location states
	full->x.sub_range(0,nL) = x;
	full->X.sub_matrix(0,nL,0,nL) = X;
	full->init();
}

void Kalman_SLAM::predict( BF::Linear_predict_model& lpred )
{
		// build full predict model by augmenting location predict with map identity predict model
	BF::Linear_predict_model fullpred (full->x.size(), lpred.q.size());
	FM::identity (fullpred.Fx);
	fullpred.Fx.sub_matrix(0,lpred.Fx.size1(),0,lpred.Fx.size2()) = lpred.Fx;
	fullpred.G.clear ();
	fullpred.G.sub_matrix(0,lpred.G.size1(), 0,lpred.G.size2()) = lpred.G;
	fullpred.q = lpred.q;
		// predict 
	full->predict (fullpred);
}

void Kalman_SLAM::observe( unsigned feature, const Feature_observe& fom, const FM::Vec& z )
{
	const size_t z_size = z.size();
		// size consistency, single state feature
	if (fom.Hx.size1() != z_size)
		error (BF::Logic_exception("observation and model size inconsistent"));
	if (fom.Hx.size2() != nL+1)
		error (BF::Logic_exception("fom and location  size inconsistent"));
		
	if (feature >= nM) {
		error (BF::Logic_exception("Observe non existing feature"));
		return;
	}
	// TODO Implement nonlinear form
		// create a augmented sparse observe model for full states
	BF::Linear_uncorrelated_observe_model fullobs(full->x.size(), z_size);
	fullobs.Hx.clear();
	fullobs.Hx.sub_matrix(0,z_size, 0,nL) = fom.Hx.sub_matrix(0,z_size, 0,nL);
	fullobs.Hx.sub_matrix(0,z_size, nL+feature,nL+feature+1) = fom.Hx.sub_matrix(0,z_size, nL,nL+1);
	fullobs.Zv = fom.Zv;
	full->observe(fullobs, z);
}

void Kalman_SLAM::observe_new( unsigned feature, const Feature_observe_inverse& fom, const FM::Vec& z )
// fom: must have a the special from required for SLAM::obeserve_new
// Feature numbers must be sequential for efficiency, otherwise state is sparse
{
		// size consistency, single state feature
	if (fom.Hx.size1() != 1)
		error (BF::Logic_exception("observation and model size inconsistent"));
		
		// make new filter with additional feature state
	if (feature >= nM)
	{
		nM = feature+1;	
		Kalman_filter_generator::Filter_type* nf = fgenerator.generate(nL+nM);
		FM::noalias(nf->x.sub_range(0,full->x.size())) = full->x;
		FM::noalias(nf->X.sub_matrix(0,full->x.size(),0,full->x.size())) = full->X;

		fgenerator.dispose(full);
		full = nf;
	}
		// build augmented location and observation
	FM::Vec sz(nL+z.size());
	sz.sub_range(0,nL) = full->x.sub_range(0,nL);
	sz.sub_range(nL,nL+z.size() )= z;

	// TODO use named references rather then explict Ha Hb
	FM::Matrix Ha (fom.Hx.sub_matrix(0,1, 0,nL) );
	FM::Matrix Hb (fom.Hx.sub_matrix(0,1, nL,nL+z.size()) );
	FM::Matrix tempHa (1,nL);
	FM::Matrix tempHb (1,sz.size());
		// feature state and variance
	full->x[nL+feature] = fom.h(sz)[0];
	full->X(nL+feature,nL+feature) = ( FM::prod_SPD(Ha,full->X.sub_matrix(0,nL, 0,nL),tempHa) +
													  FM::prod_SPD(Hb,fom.Zv,tempHb)
													 ) (0,0);
		// feature covariance with existing location and features
	//	full->X.sub_matrix(nL+feature,nL+feature+1,0,nL+nM) = prod(Ha,full->X.sub_matrix(0,nL, 0,nL+nM) );
	{	// ISSUE old uBLAS has problems assigning to symmetric proxy as above, go element by element
		const FM::Matrix cross (prod(Ha, full->X.sub_matrix(0,nL, 0,nL+nM)) );
		for (size_t i = 0; i != nL+nM-1; ++i)
			full->X(nL+feature, i) = cross(0,i);
	}
		
	full->init ();
}

void Kalman_SLAM::observe_new( unsigned feature, const FM::Float& t, const FM::Float& T )
// Feature numbers must be sequential for efficiency, otherwise state is sparse
{
		// make new filter with additional feature state
	if (feature >= nM)
	{
		nM = feature+1;
		Kalman_filter_generator::Filter_type* nf = fgenerator.generate(nL+nM);
		FM::noalias(nf->x.sub_range(0,full->x.size())) = full->x;
		FM::noalias(nf->X.sub_matrix(0,full->x.size(),0,full->x.size())) = full->X;

		fgenerator.dispose(full);
		full = nf;
	}
	else
	{
		full->x[nL+feature] = t;
		full->X(nL+feature,nL+feature) = T;
	}
	full->init ();
}

void Kalman_SLAM::forget( unsigned feature, bool must_exist )
{
	full->x[nL+feature] = 0.;
			// ISSUE old uBLAS has problems accessing the lower symmetry via a sub_matrix proxy, zero via upper symmetry
	// zero( full->X.sub_matrix(0,full->X.size1(), nL+feature,nL+feature+1) );
	zero( full->X.sub_matrix(0,nL+feature, nL+feature,nL+feature+1) );
	zero( full->X.sub_matrix(nL+feature,nL+feature+1, nL+feature,full->X.size1()) );
	full->init();
}

void Kalman_SLAM::statistics_sparse( BF::Kalman_state_filter& kstats ) const
{
	const size_t k = std::min(kstats.x.size(), full->x.size());
	kstats.x.clear(); kstats.X.clear();
	kstats.x.sub_range(0,k) = full->x.sub_range(0,k);
	kstats.X.sub_matrix(0,k, 0,k) = full->X.sub_matrix(0,k, 0,k);
}
	
void Kalman_SLAM::decorrelate( Bayesian_filter::Bayes_base::Float d )
// Reduce correlation by scaling cross-correlation terms
{
	size_t i,j;
	const size_t n = full->X.size1();
	for (i = 1; i < n; ++i)
	{
		FM::SymMatrix::Row Xi(full->X,i);
		for (j = 0; j < i; ++j)
		{
			Xi[j] *= d;
		}
		for (j = i+1; j < n; ++j)
		{
			Xi[j] *= d;
		}
	}
	full->init();
}

}//namespace SLAM

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 */

namespace SLAM_filter
{
namespace BF = Bayesian_filter;
namespace FM = Bayesian_filter_matrix;

class SLAM : public BF::Bayes_filter_base
/*
 * SLAM : Simulataneous Location and Mapping
 *  Abstract representation of general SLAM
 * The abstraction  represents the feature observation functions
 *  Observe parameters are defined:
 *   feature: A arbitary unique number to label each feature in the map.
 *   fom: feature observe model.
 * A complete SLAM solution must also represent predict; of location only or location
 * and map. This is not include in the abstraction as no single implementation can deal with
 * a general stochastic predict model.
 */
{
public:
	SLAM ()
	{}

									// Observation models
	typedef BF::Linrz_uncorrelated_observe_model Feature_observe;
	// Linearised observation model
	//  Observation z = h(lt) where lt is the vector of location state augmented with the associated feature state

	typedef BF::Uncorrelated_addative_observe_model Feature_observe_inverse;
	// Inverse model required for observe_new
	//  Feature state t = h(lz)	where lz is the vector of location state augmented with observation z

									// Single feature observe (single element vectors)
	virtual void observe( unsigned feature, const Feature_observe& fom, const FM::Vec& z ) = 0;
	// Feature observation (fuse with existing feature)
	virtual void observe_new( unsigned feature, const Feature_observe_inverse& fom, const FM::Vec& z ) = 0;
	// New Feature observation (new or overwrite existing feature)
	virtual void observe_new( unsigned feature, const FM::Vec& t, const FM::Vec& T ) = 0;
	// New Feature directly from Feature statistics: mean t and variance T (overwrite existing feature)

	virtual void forget( unsigned feature, bool must_exist = true ) = 0;
	// Forget information associated with a feature: feature number can be reused for a new feature

										// Multi feature observe: Defaults using multiple single feature observes
	typedef std::vector<unsigned> afeatures_t;	// Use a vector to store feature assocations
	virtual void multi_observe( afeatures_t& features, const Feature_observe& fom, const FM::Vec& z )
	{
		error (BF::Logic_exception("Unimplemented"));
	}
	virtual void multi_observe_new( afeatures_t& features, const Feature_observe_inverse& fom, const FM::Vec& z )
	{
		error (BF::Logic_exception("Unimplemented"));
	}
	virtual void multi_observe_new( afeatures_t& features, const FM::Vec& t, const FM::Vec& T )
	{
		error (BF::Logic_exception("Unimplemented"));
	}
	virtual void multi_forget( afeatures_t& features, bool must_exist = true )
	{
		error (BF::Logic_exception("Unimplemented"));
	}

	virtual void update () = 0;
	// Update SLAM state after a sequence of observe or forget
};


}//namespace SLAM

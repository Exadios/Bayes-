/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 *  Kalman filter representing representation of SLAM
 *  A VERY simple fixed size implementation
 *
 * Reference
 *  [1] "A Solution to the Simultaneous Localization and Map Building (SLAM) Problem"
 *   MWM Gamini Dissanayake, Paul Newman, Steven Clark, Hugh Durrant-Wyte, M Csorba, IEEE T Robotics and Automation vol.17 no.3 June 2001
 */


namespace SLAM_filter
{


class Kalman_SLAM : public SLAM
{
public:
	Kalman_SLAM( Bayesian_filter::Linrz_kalman_filter& full_filter );
	void predict( Bayesian_filter::Linear_predict_model& m ); 	// TODO Allow nonlinear

	void observe( unsigned feature, const Feature_observe& fom, const FM::Vec& z);
	void observe_new( unsigned feature, const Feature_observe_inverse& fom, const FM::Vec& z);
	void observe_new( unsigned feature, const FM::Vec& t, const FM::Vec& T);
	void forget( unsigned feature, bool must_exist = true);

	void update()
	// Compute sample mean and covariance statistics of filter
	{
		full.update();
	}
	void decorrelate( Bayesian_filter::Bayes_base::Float d);

protected:
	Bayesian_filter::Linrz_kalman_filter& full;		// Full Kalman representation of state. Reference to filter parameter in constructor

private:
	unsigned nL;		// No of location states
	unsigned nM;		// No of map states
};


}//namespace SLAM

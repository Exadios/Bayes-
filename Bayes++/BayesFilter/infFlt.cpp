/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Filter.
 */
#include "matSup.h"
#include "infFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Information_filter::Information_filter (size_t x_size, size_t z_initialsize) :
		Extended_filter(x_size),
		y(x_size), Y(x_size,x_size),
		i(x_size), I(x_size,x_size),
		ZI(Empty)
/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
	update_required = true;	// Not a valid state, init is required before update can be used
}

Information_filter& Information_filter::operator= (const Information_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	y = a.y;
	Y = a.Y;
	return *this;
}


void Information_filter::init ()
/*
 * Initialise the filter from x,X
 * Precondition:
 *		x, X
 * Postcondition:
 *		x, X is PSD
 *		y, Y is PSD
 */
{
						// Information
	Float rcond = UdUinversePD (Y, X);
	rclimit.check_PD(rcond, "Xi not PD");
						// Information state
	y.assign (prod(Y,x));
	update_required = false;
}

void Information_filter::init_information (const Vec& yi, const SymMatrix& Yi)
/*
 * Special Initialisation directly from Information yi,Yi
 * Postcondition:
 *		x, X is PSD
 *		y, Y is PSD
 */
{
	y = yi; Y = Yi;
	update_required = true;
	update();		// Enforce postcond
}

void Information_filter::update ()
/*
 * Recompute x,X from y,Y
 *  Optimised using update_required (postcondition met iff update_required false)
 * Precondition:
 *		y, Y
 * Postcondition:
 *		x=X*y, X=inv(Y) is PSD
 *		y, Y is PSD
 */
{
	if (update_required)
	{		// Covariance
		Float rcond = UdUinversePD (X, Y);
		rclimit.check_PD(rcond, "Y not PD in update");

		x.assign (prod(X,y));
		update_required = false;
	}
}

Bayes_base::Float
 Information_filter::predict (Linrz_predict_model& f)
/*
 * Extented linrz information prediction
 *  Computation is through state to accommodate linearied model
 */
{
	update ();			// x,X required
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict information matrix, and state covariance
	RowMatrix temp(f.Fx.size1(), X.size2());
	X = prod_SPD(f.Fx,X, temp) + prod_SPD(f.G, f.q);

						// Information
	Float rcond = UdUinversePD (Y, X);
	rclimit.check_PD(rcond, "X not PD in predict");
						// Predict information state
	y.assign (prod(Y,x));
	return rcond;
}

inline void Information_filter::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		ZI.resize(z_size,z_size);
	}
}

Bayes_base::Float
 Information_filter::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* correlated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

	Vec zz = s + prod(h.Hx,x);		// Strange EIF obsevation object object

						// Observation Information
	Float rcond = UdUinversePD (ZI, h.Z);
	rclimit.check_PD(rcond, "Z not PD in observe");
												// Calculate EIF i
	i = prod(trans(h.Hx), prod(ZI,zz));
												// Calculate EIF I
	RowMatrix temp = prod(ZI, h.Hx);
	I = prod(trans(h.Hx), temp );

	y += i;
	Y += I;
	update_required = true;

	assert_isPSD (Y);
	return rcond;
}

Bayes_base::Float
 Information_filter::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
/* Extended linrz uncorrelated observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

	Vec zz = s + prod(h.Hx,x);		// Strange EIF obsevation object object

						// Observation Information
	Float rcond = UdUrcond_vec(h.Zv);
	rclimit.check_PD(rcond, "Zv not PD in observe");
	ZI.clear();
	for (size_t w = 0; w < h.Zv.size(); ++w)
		ZI(w,w) = 1./ h.Zv[w];					// inverse(Z)
												// Calculate EIF i
	i = prod(trans(h.Hx), prod(ZI,zz));						// ISSUE: Efficiency ZI is diagonal
												// Calculate EIF I
	RowMatrix temp = prod(ZI, h.Hx);
	I = prod(trans(h.Hx), temp );

	y += i;
	Y += I;
	update_required = true;

	assert_isPSD (Y);
	return rcond;
}


Information_joseph_filter::Information_joseph_filter (size_t x_size, size_t z_initialsize) :
		Information_filter(x_size, z_initialsize),
		t(x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
}

Information_joseph_filter::Predict_temp::Predict_temp (size_t x_size) :
/* Construct intermediate space for predict
 */
	inv_Q(x_size, x_size),
	A(x_size, x_size),
	inv_AQ(x_size, x_size),
	Chi(x_size, x_size),
	IChi(x_size, x_size),
	Ywork(x_size)
{
}	

void Information_joseph_filter::predict (Linear_invertable_predict_model& f)
/*
 * Linear information prediction
 *  Computation is through information state only
 *  Uses x(k+1|k) = Fx * x(k|k) instead of extended x(k+1|k) = f(x(k|k))
 * Prediction is done completely on y,Y
 * Requires y(k|k), Y(k|k)
 * Predicts y(k+1|k), Y(k+1|k)
 */
{
						// Inverse Prediction noise
	t.inv_Q.clear();
	mult_SPDT (f.inv.G, f.inv.q, t.inv_Q);
	
						// Joseph Information update using, f.inv.Fx
	t.A.clear();
	mult_SPDT (f.inv.Fx, Y, t.A, t.Ywork);

	Float rcond = UdUinversePD (t.inv_AQ, t.A+t.inv_Q);
	rclimit.check_PD(rcond, "(Fx'.Y.Fx+inv_Q) not PD in predict");

	t.Chi.assign (prod(t.A,t.inv_AQ));
	FM::identity(t.IChi); t.IChi -= t.Chi;

	RowMatrix temp1(t.IChi.size1(), t.A.size2()), temp2(t.Chi.size1(), t.inv_Q.size2());
	Y = prod_SPD(t.IChi,t.A, temp1) + prod_SPD(t.Chi, t.inv_Q, temp2);
						// Information state
	y = prod(prod(t.IChi,trans(f.inv.Fx)), y);

	update_required = true;
}

}//namespace

/*
 * MODULE : bayesFilter - Bayesian Filtering
 * Provides a MuPAD module interface to a filtering class
 */

// Shut up VC++, it has its own problems!
#if defined(_MSC_VER) && _MSC_VER <= 1200
#pragma warning(disable:4786)	// indentifier more 255 chars
#endif

#include <typeinfo.h>
#include <string>
#include <map>
								// MF macro is define by MuPAD module system
#undef MF
#include "BayesFilter/allFlt.h" // Include all of Bayesian Filtering library
#include <boost/random.hpp>		// Fast and good random numbers
#include "TNear.h"				// Nearest neighbour
#include "MuPadConvert.h"

MMG( info = "Module: Bayesian Filtering" ) 

namespace MuC = MuPAD_convert;
namespace FM = Bayesian_filter_matrix;


class MuFilter
/*
 * MuPAD Runtime abstract filter
 *  Polymorphic filter created at Runtime
 *  By forceing all filter classes used at runtime we avoid the need for RTTI
 *  support in filtering library
 */
{
public:
	virtual ~MuFilter() = 0 {}; 
	void checkdef(const Bayesian_filter::Bayes_filter_base* derived) const;
};


void MuFilter::checkdef(const Bayesian_filter::Bayes_filter_base* derived) const
/*
 * Check that derived is defined. If it isn't generate MFerror based on base name
 */
{
	static std::string error;
	if (!derived) {
		error = std::string("Operation undefine on filter type: ") + typeid(*this).name();
		MFerror( error.c_str() );
	}
}

struct Mu_SIR : public MuFilter, public Bayesian_filter::SIR_filter
{
	Mu_SIR(unsigned x, unsigned s, Bayesian_filter::SIR_random& r) :
		SIR_filter(x,s,r) {}
};

struct Mu_SIR_kalman : public MuFilter, public Bayesian_filter::SIR_kalman_filter
{
	Mu_SIR_kalman(unsigned x, unsigned s, Bayesian_filter::SIR_random& r) :
		SIR_kalman_filter(x,s,r) {}
};

struct Mu_Unscented : public MuFilter, public Bayesian_filter::Unscented_filter
{
	Mu_Unscented(unsigned x) : Unscented_filter(x) {}
};

struct Mu_Cov : public MuFilter, public Bayesian_filter::Covariance_filter
{
	Mu_Cov(unsigned x) : Covariance_filter(x) {}
};

struct Mu_Inf : public MuFilter, public Bayesian_filter::Information_filter
{
	Mu_Inf(unsigned x) : Information_filter(x) {}
};

struct Mu_InfJo : public MuFilter, public Bayesian_filter::Information_joseph_filter
{
	Mu_InfJo(unsigned x) : Information_joseph_filter(x) {}
};

struct Mu_SRIF : public MuFilter, public Bayesian_filter::Information_root_filter
{
	Mu_SRIF(unsigned x) : Information_root_filter(x) {}
};

struct Mu_UD : public MuFilter, public Bayesian_filter::UD_filter
{
	Mu_UD(unsigned x) : UD_filter(x, x) {}
};



class BFilter_handler
/*
 * MuPAD interface for Bayesian Filter object
 * Allow them to be accessed by integer handles
 *  Handles start at 1 and increment with each make
 * TODO: Balance the map by using random handles
 */
{
public:
	BFilter_handler();
	int make (MuFilter*);
	MuFilter* remove (MTcell arg);
	MuFilter* get_filter(MTcell arg);
private:
	typedef std::map<int,MuFilter*> Flt_map;
	Flt_map fltmap;
	unsigned handle;
};

BFilter_handler::BFilter_handler()
{
	handle = 0;
}

int BFilter_handler::make (MuFilter* f)
/*
 * Make a handle for a filter
 * Handle overflow by warning and return 0 if handle already in use;
 */
{
	++handle;
	if (handle == 0)
	{
		MFputs( "Filter handle overflow: may prevent creation of more filters" );
		++handle;
	}

	std::pair<Flt_map::iterator, bool> i = fltmap.insert(std::make_pair(handle, f));
	if (i.second)
		return handle;
	else
	{	// Not inserted
		MFputs( "Could create unique handle for filter: no filter created" );
		return 0;
	}
};

MuFilter*  BFilter_handler::remove (MTcell arg)
/*
 * Remove a handle for a filter and return filter pointer if valid
 */
{
	if (!MFisInt(arg))
		MFerror( "Filter handle not an Int" );
	const int h = MFint(arg);
	if (h <=0) {
		MFerror( "Filter handle incorrect" );
		return 0;	// Never reached
	}
	Flt_map::iterator e = fltmap.find(h);
	if (e ==fltmap.end()) {
		MFerror( "Filter handle incorrect" );
		return 0;	// Never reached
	}

	MuFilter* f = (*e).second;
	fltmap.erase(e);

	return f;
}

MuFilter* BFilter_handler::get_filter (MTcell arg)
/*
 * Convert a handle into a filter pointer
 */
{
	if (!MFisInt(arg))
		MFerror( "Filter handle not an Int" );
	const int h = MFint(arg);
	if (h <=0) {
		MFerror( "Filter handle incorrect" );
		return 0;	// Never reached
	}
	Flt_map::iterator e = fltmap.find(h);
	if (e ==fltmap.end()) {
		MFerror( "Filter handle incorrect" );
		return 0;	// Never reached
	}

	return (*e).second;
}

// Maintain state of handles functions innovations
BFilter_handler Handles;



class Boost_random
/*
 * Random number distributions
 */
{
public:
	Boost_random() : gen_normal(rng), gen_uniform(rng)
	{
	}
	double normal(const double mean, const double sigma)
	{
		boost::normal_distribution<boost::mt19937> gen(rng, mean, sigma);
		return gen();
	}
	void normal(FM::Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_normal);
	}
	void uniform_01(FM::Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_uniform);
	}
private:
	boost::mt19937 rng;
	boost::normal_distribution<boost::mt19937> gen_normal;
	boost::uniform_01<boost::mt19937> gen_uniform;
};

// Maintain state of random numbers between function innovations
static Boost_random Random;


class MuPAD_function_model : public Bayesian_filter::Function_model
/*
 * Function model rapper for a MuPad function
 */
{
public:
	MuPAD_function_model(const MTcell ff);
	virtual const FM::Vec& fx(const FM::Vec& x) const;
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.
private:
	const MTcell Mu_fn;
	mutable FM::Vec rfx;
};

MuPAD_function_model::MuPAD_function_model(const MTcell ff) :
	Mu_fn(ff), rfx(FM::Empty)
{}

const FM::Vec& MuPAD_function_model::fx(const FM::Vec& x) const
{
	MTcell xarray = MuC::Array(x);
	MTcell fxarray = MFcall( MFcopy(Mu_fn), 1, xarray );

	if (MFisExpr(fxarray))
		MFerror ("prediction function not found");	// ISSUE: MFerror may leak arrays

	rfx = MuC::Vector(fxarray);
	MFfree(xarray);
	MFfree(fxarray);
	return rfx;
};



class Likelihood_observe_MuPAD : public Bayesian_filter::Likelihood_observe_model
/*
 * Likelihood observe model using a MuPad expression
 */
{
public:
	Likelihood_observe_MuPAD(unsigned z_size, const MTcell f);
	~Likelihood_observe_MuPAD();
	Float L(const FM::Vec& x) const;
	void Lz(const FM::Vec& z);
private:
	const MTcell Mu_fn;
	MTcell Mu_z;
};

Likelihood_observe_MuPAD::Likelihood_observe_MuPAD(unsigned z_size, const MTcell f) :
	Likelihood_observe_model(z_size),
	Mu_fn(f)
{
	Mu_z = 0;
}

Likelihood_observe_MuPAD::~Likelihood_observe_MuPAD()
{
	if (Mu_z)
		MFfree(Mu_z);
}

Likelihood_observe_MuPAD::Float Likelihood_observe_MuPAD::L(const FM::Vec& x) const
{
	MTcell Mu_x = MuC::Array(x);
	MTcell Mu_L = MFcall( MFcopy(Mu_fn), 2, MFcopy(Mu_z), Mu_x );

	if (MFisExpr(Mu_L))
		MFerror ("Likelihood observe function not found");	// ISSUE: MFerror may leak arrays

	Float L = MFdouble(Mu_L);		// ISSUE: Assume Float is double
	MFfree(Mu_x);
	MFfree(Mu_L);
	return L;
};

void Likelihood_observe_MuPAD::Lz(const FM::Vec& z)
{
	if (Mu_z)
		MFfree(Mu_z);
	Mu_z = MuC::Array(z);
};




class Boost_SIR_random_helper : public Bayesian_filter::SIR_filter::Random
/*
 * Random number generator for SIR_filter
 * Uses global Random
 */
{
public:
	Boost_SIR_random_helper()
	{}
	void normal(FM::Vec& v)
	{
		::Random.normal(v);
	}
	void uniform_01(FM::Vec& v)
	{
		::Random.uniform_01(v);
	}
};

class Boost_Predict_random_helper : public Bayesian_filter::General_LiAd_predict_model::Random
/*
 * Random number generator for General_LiAd_predict_model
 * Uses global Random
 */
{
public:
	void normal(FM::Vec& v)
	{
		::Random.normal(v);
	}
};



class Filter_maker
/*
 * Make instances of Bayes filters by name
 */
{
public:
	Filter_maker();
	typedef MuFilter* Filter;
	Filter make (MTcell args);
private:
	MTcell MVargs;		// Same name as in MuPAD function so macros work!
	unsigned state_size();
	Boost_SIR_random_helper randomHelper;

	typedef Filter (Filter_maker::* pMake_mem)();
	typedef std::map<std::string, pMake_mem> Make_map;
	Make_map filterMakers;

	// Filter makers
	Filter make_SIR();
	Filter make_SIR_Kalman();
	Filter make_Unscented();
	Filter make_Covariance();
	Filter make_Information();
	Filter make_Information_joseph();
	Filter make_SRIF();
	Filter make_UD();
};


Filter_maker::Filter_maker()
/*
 * Construct name-> make mapping
 */
{
	filterMakers["SIR"] = &Filter_maker::make_SIR;
	filterMakers["SIR_Kalman"] = &Filter_maker::make_SIR_Kalman;
	filterMakers["Unscented"] = &Filter_maker::make_Unscented;
	filterMakers["Covariance"] = &Filter_maker::make_Covariance;
	filterMakers["Information"] = &Filter_maker::make_Information;
	filterMakers["InformationJoseph"] = &Filter_maker::make_Information_joseph;
	filterMakers["SRIF"] = &Filter_maker::make_SRIF;
	filterMakers["UD"] = &Filter_maker::make_UD;
}

Filter_maker::Filter
 Filter_maker::make (MTcell args)
/*
 * Use arg to Make instance of a filter
 */
{
	MVargs = args;					// Allow MFarg macros to work in instance
	if (MVnargs < 1)
		MFerror("No filter name argument");
	MFargCheck (1, DOM_STRING);

	std::string filterTypeName(MFstring(MFarg(1)));
								// Find the entry in filterMakers map
	Make_map::iterator i = filterMakers.find (filterTypeName);

	Filter f = 0;
	if (i != filterMakers.end())
	{	// Call the make member function found to create a filter
		pMake_mem pmake = i->second;
		f = (this->*pmake)();
	}
	return f;
};
	


unsigned Filter_maker::state_size()
{
	MFargCheck (2, DOM_INT);
	int ss = MFint(MFarg(2));
	if (!(ss > 0))
		MFerror ( "state_size must > 0" );
	return unsigned(ss);
}

Filter_maker::Filter
 Filter_maker::make_SIR()
{
	MFnargsCheckRange (4,6);

				// Initialise Sample S with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));
				// Defaultable parameters
	int sample_size = 1000;
	double rougheningK = -1;
	if (MVnargs >=5)
	{
		MFargCheck (5, DOM_INT);
		sample_size = MFint(MFarg(5));
		if (!(sample_size > 0))
			MFerror ( "sample_size must > 0" );
	}
	if (MVnargs >=6)
	{
		MFargCheck (6, DOM_FLOAT);
		rougheningK = MFdouble(MFarg(6));
		if (!(rougheningK >= 0))
			MFerror ( "rougheningK must >= 0" );
	}

	Mu_SIR* f = new Mu_SIR(state_size(), sample_size, randomHelper);

	// Use a temporary filter to create samples
	Bayesian_filter::SIR_kalman_filter tempS(f->S.size1(), f->S.size2(), randomHelper);
	tempS.init_kalman (init_x,init_X);
	f->init(tempS.S);
				// Overide filter default
	if (rougheningK != -1.)
		f->rougheningK = rougheningK;

	return f;
}

Filter_maker::Filter
 Filter_maker::make_SIR_Kalman()
{
	MFnargsCheckRange (4,6);

				// Initialise Sample S with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));
				// Defaultable parameters
	int sample_size = 1000;
	double rougheningK = -1;
	if (MVnargs >=5)
	{
		MFargCheck (5, DOM_INT);
		sample_size = MFint(MFarg(5));
		if (!(sample_size > 0))
			MFerror ( "sample_size must > 0" );
	}
	if (MVnargs >=6)
	{
		MFargCheck (6, DOM_FLOAT);
		rougheningK = MFdouble(MFarg(6));
		if (!(rougheningK >= 0))
			MFerror ( "rougheningK must >= 0" );
	}

	Mu_SIR_kalman* f = new Mu_SIR_kalman(state_size(), sample_size, randomHelper);

	f->init_kalman (init_x,init_X);
				// Overide filter default
	if (rougheningK != -1.)
		f->rougheningK = rougheningK;

	return f;
}

Filter_maker::Filter
 Filter_maker::make_Unscented()
{
	MFnargsCheck (4);

	// Initialise with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));
	
	Mu_Unscented* f = new Mu_Unscented(state_size());
	f->init_kalman (init_x,init_X);
	return f;
}

Filter_maker::Filter
 Filter_maker::make_Covariance()
{
	MFnargsCheck (4);

	// Initialise with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));

	Mu_Cov* f = new Mu_Cov(state_size());
	f->init_kalman (init_x,init_X);
	return f;
}

Filter_maker::Filter
 Filter_maker::make_Information()
{
	MFnargsCheck (4);

	// Initialise with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));

	Mu_Inf* f = new Mu_Inf(state_size());
	f->init_kalman (init_x,init_X);
	return f;
}

Filter_maker::Filter
 Filter_maker::make_Information_joseph()
{
	MFnargsCheck (4);

	// Initialise with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));

	Mu_InfJo* f = new Mu_InfJo(state_size());
	f->init_kalman (init_x,init_X);
	return f;
}

Filter_maker::Filter
 Filter_maker::make_SRIF()
{
	MFnargsCheck (4);

	// Initialise with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));

	Mu_SRIF* f = new Mu_SRIF(state_size());
	f->init_kalman (init_x,init_X);
	return f;
}

Filter_maker::Filter
 Filter_maker::make_UD()
{
	MFnargsCheck (4);

	// Initialise Sample S with state and covariance
	FM::Vec init_x = MuC::Vector(MFarg(3));
	FM::Matrix init_X = MuC::Matrix(MFarg(4));

	Mu_UD* f = new Mu_UD(state_size());
	f->init_kalman (init_x,init_X);
	return f;
}


// Maintain a set of default filters to make
Filter_maker DefaultFilters;


//
// Module functions definitions
//

MFUNC( bfilter, MCnop )
/*
 * DOM_INT bfilter(String filterTypeName, DOM_INT nstate, DOM_INT nsamples)
 * Create a filter from its names type and return its handle (Int)
 */
{
	// Construct Filter  - Must catch exceptions for MuPAD kernel
	MuFilter* filter;
	try
	{
		filter = DefaultFilters.make (MVargs);
		if (!filter)
			MFerror( "Unknown filter type name" );
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
	catch (...)
	{
		MFerror( "Caused a system exception" );
	}

	// Make a handle to return.
	const int h = Handles.make (filter);
	if (!h)		// Not Created
		delete filter;
	MFreturn (MFint(h));
} MFEND 


MFUNC( free, MCnop )
/*
 * DOM_NULL free(DOM_INT hFtiler)
 * free the resources associated with a filter handle
 * The handle thus invalidated
 */
{
	MFnargsCheck( 1 );
	MuFilter* filter = Handles.remove(MFarg(1));

	if (filter)
	{
		delete filter;
	}

	MFreturn (MFcopy(MVnull));
} MFEND


MFUNC( predict_functional, MCnop )
/*
 * DOM_NULL predict_functional (DOM_INT hFilter, DOM_PROC f)
 *  Predict filter with zero noise through function f(x)
 */
{
	MFnargsCheck (2);
	MuFilter* filter = Handles.get_filter(MFarg(1));

	// Must be a Functional filter
	Bayesian_filter::Functional_filter* f = dynamic_cast<Bayesian_filter::Functional_filter*>(filter);
	filter->checkdef(f);

	class FF : public Bayesian_filter::Functional_predict_model
	{	// Functional form
		Bayesian_filter::Function_model& ff;
	public:
		FF (Bayesian_filter::Function_model& f_init) :
			ff(f_init)
		{}
		virtual const FM::Vec& fx(const FM::Vec& x) const
		{	return ff.fx(x);
		}
	};

	try
	{
		// Predict filter using MuPAD expression
		MuPAD_function_model mu_f(MFarg(2));

		f->predict (FF(mu_f));
		MFreturn (MFcopy(MVnull));
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND


MFUNC( predict_additive, MCnop )
/*
 * DOM_REAL predict_additive (DOM_INT hFilter, DOM_PROC f, DOM_Array q, DOM_Array G)
 *  Predict filter with additive noise Gq
 */
{
	MFnargsCheck (4);
	MuFilter* filter = Handles.get_filter(MFarg(1));
	FM::Vec    q = MuC::Vector(MFarg(3));
	FM::Matrix G = MuC::Matrix(MFarg(4));

	// Check matrix conformance
	{
		Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
		filter->checkdef(f);

		if (f->x.size() != G.size1())
			MFerror("Mismatch in x and G size");
		if (q.size() != G.size2())
			MFerror("Mismatch in q and G size");
	}

	bool bFilterOp = false;			// Flag operation complete to avoid MFerror in try block
	double rcond;
	try
	{
		// Predict filter using MuPAD expression
		MuPAD_function_model mu_f(MFarg(2));
		Bayesian_filter::Simple_additive_predict_model model(mu_f, G, q);
		MFerror ("Not implemented");

#ifdef UNIMPLEMENTED
// TODO Requires General_LzAd_predict_model
		Boost_Predict_random_helper predict_random;
		Bayesian_filter::General_LzAd_predict_model Lmodel(G.size1(), q.size(), predict_random);
		// Use mu_f to construct Lmodel
		FM::copy (model.q, Lmodel.q);
		FM::copy (model.G, Lmodel.G);

		// Must be a Sample filter
		if (Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter)) {
			bFilterOp = true;
			f->predict (Lmodel);
			rcond = 0.;
		}
#endif
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}

	if (bFilterOp) {
		MFreturn (MFdouble(rcond));
	}
	else {
		filter->checkdef(0);
	}
} MFEND


MFUNC( predict_linear, MCnop )
/*
 * DOM_REAL predict_linear (DOM_INT hFilter, DOM_Array F, DOM_Array G, DOM_Array q)
 *  Predict filter with additive noise Gq
 */
{
	MFnargsCheck (4);
	MuFilter* filter = Handles.get_filter(MFarg(1));
	FM::Matrix Fx = MuC::Matrix(MFarg(2));
	FM::Matrix G = MuC::Matrix(MFarg(3));
	FM::Vec    q = MuC::Vector(MFarg(4));

	// Check matrix conformance
	{
		Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
		filter->checkdef(f);

		if (f->x.size() != Fx.size1())
			MFerror("Mismatch in x and Fx matrix");
		if (Fx.size1() != Fx.size2())
			MFerror("Mismatch Fx not square");
		if (f->x.size() != G.size1())
			MFerror("Mismatch in x and G size");
		if (q.size() != G.size2())
			MFerror("Mismatch in q and G size");
	}

	bool bFilterOp = false;			// Flag operation complete to avoid MFerror in try block
	double rcond;
	try
	{
		// Predict filter 
		Boost_Predict_random_helper random;
		Bayesian_filter::General_LiAd_predict_model model(Fx.size1(), q.size(), random);
		model.Fx = Fx;
		model.q = q;
		model.G = G;

		// Must be a Likelihood, or Sample filter
		if (Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter)) {
			bFilterOp = true;
			f->predict (model);
			rcond = 0.;
		}
		else
		if (Bayesian_filter::Linrz_filter* f = dynamic_cast<Bayesian_filter::Linrz_filter*>(filter)) {
			bFilterOp = true;
			rcond = f->predict (model);
		}
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}

	if (bFilterOp) {
		MFreturn (MFdouble(rcond));
	}
	else {
		filter->checkdef(0);
	}
} MFEND


MFUNC( observe_likelihood, MCnop )
/*
 * DOM_REAL observe_likelihood (DOM_INT hFilter, DOM_PROC L, DOM_ARRAY z)
 *  Observe filter with Likelihood L(z,x)
 */
{
	MFnargsCheck (3);
	MuFilter* filter = Handles.get_filter(MFarg(1));
	FM::Vec    z = MuC::Vector(MFarg(3));

	// Must be a Sample filter
	Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
	filter->checkdef(f);

	try
	{
		// Observe filter using MuPAD expression
		Likelihood_observe_MuPAD model(z.size(), MFarg(2));
		f->observe(model, z);
		double rcond;
		f->update (rcond);
		MFreturn (MFdouble(rcond));
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND


MFUNC( observe_linear_uncorrelated, MCnop )
/*
 * DOM_REAL observe_linear_uncorrelated (DOM_INT hFilter, DOM_ARRAY H, DOM_ARRAY z, DOM_ARRAY Zd)
 *  Observe z with noise Zv through linear H
 */
{
	MFnargsCheck (4);
	MuFilter* filter = Handles.get_filter(MFarg(1));
	FM::Matrix H = MuC::Matrix(MFarg(2));
	FM::Vec    z = MuC::Vector(MFarg(3));
	FM::Vec	  Zv = MuC::Vector(MFarg(4));

	{
		const Bayesian_filter::State_filter* f = dynamic_cast<const Bayesian_filter::State_filter*>(filter);
		filter->checkdef(f);

		// Check matrix conformance
		if (f->x.size() != H.size2())
			MFerror("Mismatch in x and H matrix");
		if (z.size() != H.size1())
			MFerror("Mismatch in z and H matrix");
		if (Zv.size() != z.size())
			MFerror("Mismatch in z and Zv size");
	}

	bool bFilterOp = false;			// Flag operation to avoid MFerror in try block
	double rcond;
	try
	{
		// Observe filter
		Bayesian_filter::General_LiUnAd_observe_model model(H.size2(), z.size());
		model.Hx = H;
		model.Zv = Zv;

		// Must be a Sample, Linrz filter
		if (Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter)) {
			bFilterOp = true;
			f->observe (model, z);
			double rcond;
			f->update (rcond);
		}
		else
		if (Bayesian_filter::Linrz_filter* f = dynamic_cast<Bayesian_filter::Linrz_filter*>(filter)) {
			bFilterOp = true;
			rcond = f->observe (model, z);
		}
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}

	if (bFilterOp) {
		MFreturn (MFdouble(rcond));
	}
	else {
		filter->checkdef(0);
	}
} MFEND


MFUNC( observe_linear_correlated, MCnop )
/*
 * DOM_REAL observe_linear_correlated (DOM_INT hFilter, DOM_ARRAY H, DOM_ARRAY z, DOM_ARRAY Z)
 *  Observe z with noise Zv through linear H
 */
{
	MFnargsCheck (4);
	MuFilter* filter = Handles.get_filter(MFarg(1));
	FM::Matrix H = MuC::Matrix(MFarg(2));
	FM::Vec    z = MuC::Vector(MFarg(3));
	FM::Matrix Z = MuC::Matrix(MFarg(4));

	{
		const Bayesian_filter::State_filter* f = dynamic_cast<const Bayesian_filter::State_filter*>(filter);
		filter->checkdef(f);

		// Check matrix conformance
		if (f->x.size() != H.size2())
			MFerror("Mismatch in x and H matrix");
		if (z.size() != H.size1())
			MFerror("Mismatch in z and H matrix");
		if (Z.size1() != z.size())
			MFerror("Mismatch in z and Z size");
		if (Z.size1() != Z.size2())
			MFerror("Mismatch Z not square");
	}

	bool bFilterOp = false;			// Flag operation to avoid MFerror in try block
	double rcond;
	try
	{
		// Observe filter 
		Bayesian_filter::General_LiCoAd_observe_model model(H.size2(), z.size());
		model.Hx = H;
		model.Z = Z;
		// Must be a Sample, Linrz filter
		if (Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter)) {
			bFilterOp = true;
			f->observe (model, z);
			f->update (rcond);
		}
		else
		if (Bayesian_filter::Linrz_filter* f = dynamic_cast<Bayesian_filter::Linrz_filter*>(filter)) {
			bFilterOp = true;
			rcond = f->observe (model, z);
		}
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}

	if (bFilterOp) {
		MFreturn (MFdouble(rcond));
	}
	else {
		filter->checkdef(0);
	}
} MFEND


MFUNC( mean, MCnop )
/*
 * DOM_ARRAY mean(DOM_INT hFilter)
 *  mean of filter state
 */
{
	MFnargsCheck( 1 );
	MuFilter* filter = Handles.get_filter(MFarg(1));

	// Must be a State filter
	Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
	filter->checkdef(f);

	// Process Filter  - Must catch exceptions for MuPAD kernel
	try
	{
		f->update();
		MFreturn(MuC::Array(f->x));
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND 


MFUNC( covariance, MCnop )
/*
 * DOM_ARRAY covariance (DOM_INT hFilter)
 *  covariance of filter state
 */
{
	MFnargsCheck( 1 );
	MuFilter* filter = Handles.get_filter(MFarg(1));

	// Must be a Kalman filter
	Bayesian_filter::Kalman_filter* f = dynamic_cast<Bayesian_filter::Kalman_filter*>(filter);
	filter->checkdef(f);

	// Process Filter  - Must catch exceptions for MuPAD kernel
	try
	{
		f->update();
		MFreturn(MuC::Array(f->X));
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND 


MFUNC( unique_samples, MCnop )
/*
 * DOM_INT unique_samples (DOM_INT hFilter)
 *  no of unique (different value) samples in Filter
 */
{
	MFnargsCheck( 1 );
	MuFilter* filter = Handles.get_filter(MFarg(1));

	// Must be a Sample filter
	Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
	filter->checkdef(f);

	// Process Filter  - Must catch exceptions for MuPAD kernel
	try
	{
		unsigned nsample = f->unique_samples();
		MFreturn(MFint(nsample));
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND 


MFUNC( stochastic_samples, MCnop )
/*
 * DOM_INT stochastic_samples (DOM_INT hFilter)
 *  no of unique (stochastic history) samples in Filter
 */
{
	MFnargsCheck( 1 );
	MuFilter* filter = Handles.get_filter(MFarg(1));

	// Must be a SIR filter
	Bayesian_filter::SIR_filter* f = dynamic_cast<Bayesian_filter::SIR_filter*>(filter);
	filter->checkdef(f);

	// Process Filter  - Must catch exceptions for MuPAD kernel
	try
	{
		unsigned nsample = f->stochastic_samples;
		MFreturn(MFint(nsample));
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND 


MFUNC( sample, MCnop )
/*
 * DOM_LIST sample (DOM_INT hFilter)
 *  Return samples (list of lists) reresenting filters state probability distribution.
 *  If the filter is represented by samples these direcly used
 *  otherwise a set of samples are generated to represent it
 */
{
	const unsigned GenerateSamples = 1000;
	MFnargsCheck( 1 );
	MuFilter* filter = Handles.get_filter(MFarg(1));

	// Process filter by type - Must catch exceptions for MuPAD kernel
	try
	{
		// Sample filter is easy
		{
			Bayesian_filter::Sample_filter* sf = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
			if (sf) {
				MFreturn(MuC::ListTranspose(sf->S));
			}
		}

		// Kalman filter requres more work. Create Samples from a mean and covariance of Kalman filter
		{
			Bayesian_filter::Kalman_filter* kf = dynamic_cast<Bayesian_filter::Kalman_filter*>(filter);
			if (kf) {
				kf->update();
				Boost_SIR_random_helper randomHelper;
				Bayesian_filter::SIR_kalman_filter tempS(kf->x.size(), GenerateSamples, randomHelper);
				tempS.init_kalman (kf->x, kf->X);
				MFreturn(MuC::ListTranspose(tempS.S));
			}
		}

		// Fall through if type unknown
		filter->checkdef(0);
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND 


#error This is wrong V must be a pointer like object
class V  
{	// Vector proxy and difference for Neariest Neighbour
public:
	typedef FM::ColMatrix::const_Column NNVec;

	double distance (const V& o) const
	{
		double d = 0.;
		NNVec::const_iterator oi = (o.dp).begin();
		for (NNVec::const_iterator vi = (dp).begin(); vi != (dp).end(); ++vi, ++oi) {
			const double diff = *vi - *oi;
			d += diff*diff;
		}
		return sqrt(d);
	}
	V()	// Empty for return value
	{
	}
	V( const FM::ColMatrix& colmat, std::size_t i) : dp(colmat,i)
	{
	}
	NNVec dp;	// Vec data pointer
};

FM::Vec find_nearest(const FM::ColMatrix& S, double radius, const FM::Vec& loc)
/*
 * Find the nearest neighbour to loc in S
 */
{
	// Construct NN tree of samples
	CNearTree<V> sampleTree;
	for (std::size_t si = 0; si != S.size2(); ++si) {
		sampleTree.m_fnInsert( V(S,si) );
	}
	// Create a location as same type as a sample
	FM::ColMatrix cmloc(S.size1(),1);
	FM::column(cmloc,0) = loc;

	// Find NN
	V nn;
	if (sampleTree.m_bfnNearestNeighbor(radius, nn, V(cmloc,0)) )
	{
		FM::Vec rnn(S.size1());
		rnn = nn.dp;
		return rnn;
	}
	else
	{	// No NN
		FM::Vec rnn(FM::Empty);
		return rnn;

	}
}

MFUNC( nearest_sample, MCnop )
/*
 * Find nearest neighbour in a sample
 */
{
	MFnargsCheck( 3 );
	MuFilter* filter = Handles.get_filter(MFarg(1));
	double	  radius = MFdouble(MFarg(2));
	FM::Vec   loc = MuC::Vector(MFarg(3));

	// Process filter by type - Must catch exceptions for MuPAD kernel
	try
	{
		// Sample filter is easy
		{
			Bayesian_filter::Sample_filter* sf = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
			if (sf) {
				MFreturn( MuC::List(find_nearest(sf->S, radius, loc)) );
			}
		}

		// Fall through if type unknown
		filter->checkdef(0);
	}
	catch (Bayesian_filter::Bayes_filter_exception fe)
	{
		MFerror( fe.what() );
	}
} MFEND


MFUNC( _test, MCnop )
/*
 * testing
 */
{
	MFnargsCheck( 0 );

	MFreturn (MFcopy(MVnull));

	//FM::Vec v = MuC::Vector(MFarg(2)); 
	//MFreturn (MuC::List(v));

} MFEND

#
# Makefile: Matlab MEX interface for Bayes++
#  Build with Debug compilation but release libs
#
# NOTE: VC7 Require Users mexopts.bat to increase heap using -Zm150
BAYESFILTER = ..
CLAPACKLIB = ../../CLAPACK3/Release/clapack7.lib

MEXOPTS = -I$(BAYESFILTER) -I$(BAYESFILTER)/uBLASmatSup
LIBS = $(BAYESFILTER)/VClib/BayesFilter_UBLAS_MD7.lib $(CLAPACKLIB)

all : bf_all slam_all

clean:
	del *.dll *.ilk *.pdb


bf_all : \
	bf_scheme.dll \
	bf_predict_functional.dll \
	bf_predict_addative.dll \
	bf_predict_linear.dll \
	bf_predict_linrz.dll \
	bf_observe_likelihood.dll \
	bf_observe_linear_uncorrelated.dll \
	bf_observe_linear_correlated.dll \
	bf_observe_linrz_uncorrelated.dll \
	bf_observe_linrz_correlated.dll \
	bf_mean.dll \
	bf_cov.dll  \
	bf_unique_samples.dll \
	bf_stochastic_samples.dll \
	bf_sample.dll

slam_all : \
	slam_scheme.dll slam_observe.dll
# To be added 
# slam_observe_new.dll
	
bf_scheme.dll : bf_scheme.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_scheme.cpp $(LIBS)

bf_predict_functional.dll : bf_predict_functional.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_predict_functional.cpp $(LIBS)

bf_predict_addative.dll : bf_predict_addative.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_predict_addative.cpp $(LIBS)

bf_predict_linear.dll : bf_predict_linear.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_predict_linear.cpp $(LIBS)

bf_predict_linrz.dll : bf_predict_linrz.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_predict_linrz.cpp $(LIBS)

bf_observe_likelihood.dll : bf_observe_likelihood.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_observe_likelihood.cpp $(LIBS)

bf_observe_linear_uncorrelated.dll : bf_observe_linear_uncorrelated.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_observe_linear_uncorrelated.cpp $(LIBS)

bf_observe_linear_correlated.dll : bf_observe_linear_correlated.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_observe_linear_correlated.cpp $(LIBS)

bf_observe_linrz_uncorrelated.dll : bf_observe_linrz_uncorrelated.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_observe_linrz_uncorrelated.cpp $(LIBS)

bf_observe_linrz_correlated.dll : bf_observe_linrz_correlated.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_observe_linrz_correlated.cpp $(LIBS)

bf_mean.dll : bf_mean.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_mean.cpp $(LIBS)

bf_cov.dll : bf_cov.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_cov.cpp $(LIBS)

bf_unique_samples.dll : bf_unique_samples.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_unique_samples.cpp $(LIBS)

bf_stochastic_samples.dll : bf_stochastic_samples.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_stochastic_samples.cpp $(LIBS)

bf_sample.dll : bf_sample.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) bf_sample.cpp $(LIBS)


slam_scheme.dll : slam_scheme.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) slam_scheme.cpp ../SLAM/fastSLAM.cpp ../SLAM/kalmanSLAM.cpp $(LIBS)

slam_observe.dll : slam_observe.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) slam_observe.cpp ../SLAM/fastSLAM.cpp ../SLAM/kalmanSLAM.cpp $(LIBS)

slam_observe_new.dll : slam_observe_new.cpp matlabConvert.h matlabBfilter.h
	mex $(MEXOPTS) slam_observe_new.cpp ../SLAM/fastSLAM.cpp ../SLAM/kalmanSLAM.cpp $(LIBS)

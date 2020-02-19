
## emfit toolbox 

A set of Matlab scripts which performs fixed-effects empirical Bayes batch EM
inference by first fitting priors for a set of models from each model class,
plotting the inferred posterior parameters, performing approximate model
comparison, generating surrogate data and performing some visual comparisons of
the true and surrogate data. 

A general overview over this modelling approach can be found here: 
Huys (2017): Bayesian approaches to Learning and Decision Making. 
In: Computational Psychiatry: Mathematical Modelling of Mental Disorders. Editors: Alan Anticevic and John Murray
Preprint: https://quentinhuys.com/abstr_Huys17-BayesianModellingDecisionMaking.html

Currently the following models are implemented: 

'mBasicRescorlaWagner';			% basic Rescorla-Wagner example 
'mAffectiveGoNogo';			% Guitart et al. 2012, Neuroimage 62(1):154-66 
'mProbabilisticReward';			% Huys et al., 2013 Biology of Mood & Anxiety Disorders 3:12 
'mTwostep';				% Daw et al., 2011 Neuron 69:1204
'mEffortCollins';	 		% Gold et al., 2013 Biological Psychiatry 74:130
'mPruning'; 				% Lally et al., 2017 J. Neurosci 37(42):10215
'mEffortDDM';                           % Berwian et al., 2020 JAMA Psychiatry doi:10.1001/jamapsychiatry.2019.4971

# try it out 

You can try these out by simply changing into the emfit directory within matlab
and then running

	batchRunEMfit('mAffectiveGoNogo')

which will generate an example dataset from the model, and then run a series of
models on the data, do model comparison and generate some plots. 

# fit some data 

If your data is already in the correct format (see dataformat.txt files in each
of the model folders), then you can fit the data by 

	Data = load('mydata.mat');
	batchRunEMfit('mAffectiveGoNogo',Data);

The results will be saved to a folder named fitResults in a separate .mat file
for each model. This contains
	* individual posterior parameter means 'E'
	* individual posterior parameter variances 'V'
	* group mean prior 'alpha'
	* group prior precision 'nui'
Additional aspects such as ML parameter estimates, diagnostics etc. are
contained within 'stats'. 

# add your own

Each set of models includes 
* dataformat.txt - this describes the format in which the data needs to be 
* modelList.m - a file with brief descriptions of each of the models included 
* ll***.m - the actual model likelihood files. In general, this takes the form 
  [l,dl,dsurr] = yourlikelihood(x,D,mu,nui,doprior,options);
  where 
  	* l is either the total log likelihood of the data, or the (unnormalized)
	  log posterior if doprior=1. 
	* dl is a vector of gradients of the total log likelihood with respect to each of the parameters
	* dsurr is surrogate data drawn from the model (if included)
	* x are the parameters (NB: unconstrained)
	* mu is the mean vector of the Gaussian group prior 
	* nui is the inverse covariance of the Gaussian group prior 
	* doprior is either 0 or 1. When set to 1, l and dl refer to the log
	  posterior 
	* options.generatesurrogatedata sets whether dsurr is generated or not
* generateExampleDataset.m - this generates an example dataset when called by
  batchGenerateSurrogateData.m (which in turn is called by the above
  batchRunEMfit.m)
* surrogateDataPlots.m - this file will be run by batchRunEMfit.m and should
  compare surrogate data drawn from the fitted models to the real data. 

To add your own, just add these sets of files. Easiest is probably to
copy/adapt some of the existing ones. 

# Important details 

The core functionality is contained within lib/emfit.m, which can be run
directly without the various wrapper scripts. It allows GLM-style regressions at
the group level, for instance. 

Parameters are assumed to be drawn from Gaussian priors, and the means and
variances of these are empirically estimated. Parameters are constrained to be
positive via exp(.), and within bounds via a sigmoid transform. 


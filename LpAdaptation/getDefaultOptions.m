function defopts = getDefaultOptions(N)

% Options defaults:

% maximal number of function evaluations';
defopts.MaxEval         = 1e4*(N);

% norm of lp Ball used for finding design
% center and get volume approximation
% pn: what p norm should be used  p>0
% 1 = manhatten norm
% 2 = euclidean norm -> hyperellipsoid
% >100 infinity norm -> hyperrectangle
% 0<p<1 -> starshaped, p should't be too small,  > 0.4
defopts.pn              = 2;        

% give additional arguments for oracle
defopts.oracleInopts    = [];       

% number of oracle output, first needs to be 0/1
% (unfeasible/feasible) the rest can be some values
% (e.g. fitness values) that one wants to save
defopts.nOut            = 1;        

% Plot progress while running';
defopts.Plotting        = 'on';     

% >=0, command line messages after every i-th iteration'; MaxEval mod VerboseModulo must equal 0!!!
% when 'adapt'= 'CMA', command line messages after every
% VerboseModulo/PopSize Generation
defopts.VerboseModulo   = 1e3;

% >=0, saving after every i-th iteration';MaxEval mod SavingModulo must equal 0!!!
% when 'adapt'= 'CMA', command line messages after every
% VerboseModulo/PopSize Generation
defopts.SavingModulo    = 1e2;

% [on|off] save data to file';
defopts.bSaving         = 'on';    

% save covariance matrices' ('1' results in huge files);
defopts.bSaveCov        = 1;       

% % if interested in averaged covariance matrix; average of last opts.averageCovNum covariances
% defopts.getAverageCov   = 1;       
% 
% % 1: get average over Time for radius, Volume, mean, Cov
% defopts.averageOverTime.on = 0;    
% 
% defopts.averageOverTime.meanSize = 'min(18/valP,maxMeanSize)';

% 1: save all numLast r,mu,Q,P_emp
defopts.LastSaveAll     = 0;       

% 1: save all unfeasible solutions as well
defopts.unfeasibleSave  = 0;       

% how many of numLast elements are used to get average mu and r
% should be after burn-in period
defopts.numLast         = 'max((opts.MaxEval - 1e3*N),opts.MaxEval * 0.7)'; 

% how many covariances should be used to get average covariance
defopts.averageCovNum   = 100; 

% Default options for algorithmic parameters:
% --------------------------------------------
% Hitting probability 
defopts.valP            = 1/exp(1); %1/5;

% upper bound on inverval's size over which averaging happens
defopts.maxMeanSize = 2000;
% size of moving window to get empirical hitting probability (in number of evaluations)
defopts.windowSizeEval  = 'min(110/valP,maxMeanSize)';%300; 

% Step size of the initial covariance
defopts.r               = 1;                

% Initial Cholesky matrix   %%% only one should
defopts.initQ           = eye(N);     
% Initial Covariance matrix %%% be specified
defopts.initC           = eye(N);     

% maximal allowed step size
defopts.MaxR            = Inf;              
% minimal allowed step size
defopts.MinR            = 0;                

% maximal allowed condition
defopts.MaxCond         = 1e20*N;           

% Mean adaptation weight '1/valP*N';
defopts.N_mu            = exp(1)*N;         

% Matrix adaptation weight
% equals 1/ccov1, mueff=1, because only one candidate solution 2 / ((N+1.3)^2+mueff)
defopts.N_C             = ((N+1.3)^2+1)/2;  

% Step size increase/decrease factor
defopts.beta = '3*0.2/((N+1.3)^2+valP*PopSize)';

% Expansion upon success
defopts.ss              = '1 + beta*(1-valP)'; 
% Contraction otherwise
defopts.sf              = '1 - beta*(valP)';   

% learning rate for rank-one update;
defopts.CMA.ccov1 = '3*0.2/((N+1.3)^2+mueff)';
% learning rate for rank-mu update;
defopts.CMA.ccovmu      = 'min(1-ccov1, 3*0.2*(mueff-2+1/mueff) / ((N+2)^2+mueff*0.2))'; 

%number of candidate solutions in one iteration, population size
defopts.CMA.PopSize     =  'max(4+floor(3*log(N)),floor(2/valP))';
%learning constant for the evolution path,
% which should be usually between 1/sqrt(N) and 2/(N+1). for larger cp the
% effect of the evolution path will attenuate. the backward time horizon for
% the evolution path is roughly cp^-1
defopts.CMA.cp          = 1/sqrt(N);

% 1: covariance is adapted, 0: covariance is fixed
defopts.Cadapt          = 1;        
% 1: adapt mean
defopts.madapt          = 1;        

% 1: adapt hittin probability (to get a more accurate volume estimation or to get a better design center)
% of interest if hitP_adapt == 1
defopts.hitP_adapt      = 0;        

%decreasing for volume estimation, increasing for design center
defopts.para_hitP_adapt.PVec = [0.35,0.15,0.06,0.03,0.01]; 

% 1: fixed schedule, 0: no fixed schedule, then it will be decided on the fly when to change the hitting probability
% of interest if hitP_adapt ==1 and fixed schedule
defopts.para_hitP_adapt.fixedSchedule = 1; 

%% interesting for fixed schedule (opts.para_hitP_adapt.fixedSchedule ==1)
% proportion of maxEval
% that are used for every run with a different hitting probability,
% vector should have same length as PVec, sum of this vector should ideally
% equal 1 (otherwise function evaluations are wasted)
defopts.para_hitP_adapt.maxEvalSchedule = [1/2,1/8,1/8,1/8,1/8];%[3/5,3/50,4/50,6/50,7/50]; 

% over how many
% samples of each run should be averaged to get radius r and mean mu
% of found feasible region
% of interest if hitP_adapt ==1 and no fixed schedule
defopts.para_hitP_adapt.numLastSchedule = [1/2,3/4,3/4,3/4,3/4]; 

%% interesting for variable schedule ( opts.para_hitP_adapt.fixedSchedule ==0 )
% every which iteration it should be tested if step size is in steady state
defopts.para_hitP_adapt.testEvery = 'min(18/valP,maxMeanSize)'; 

% use same testStart for all properties that are checked before hitP is adapted
% after how many iterations it should be tested for fist time; needs to be > 2*opts.para_hitP_adapt.meanSize
defopts.para_hitP_adapt.testStart = 'max([2*opts.para_hitP_adapt.stepSize.meanSize,2*opts.para_hitP_adapt.hitP.meanSize,2*opts.para_hitP_adapt.VolApprox.meanSize])' ; 

%compare two means with each other. meanSize defines the size of how how many r values are used for one mean
defopts.para_hitP_adapt.stepSize.meanSize = 'min(18/valP,maxMeanSize)';%'18/valP';%50; 
%if abs(mean1-mean2)/(mean1/2+mean2/2) < deviation, start step size adaptation
defopts.para_hitP_adapt.stepSize.deviation = 0.001; 

%loewner, axis-alligned bounding-box volume approximation
defopts.para_hitP_adapt.VolApprox.meanSize = 'min(18/valP,maxMeanSize)';
defopts.para_hitP_adapt.VolApprox.deviation = 0.001;

% hitting probability
defopts.para_hitP_adapt.hitP.meanSize = 'min(30/valP,maxMeanSize)';
defopts.para_hitP_adapt.hitP.deviation = 0.001;

% defining how many samples of each run are used for calculating the
% average if no fixed schedule for the changing of the hitting probability
%number between 0 and 1;
defopts.para_hitP_adapt.meanOfLast = 1/4;%1/2; 

defopts.para_hitP_adapt.deviation_stop = 0.01;

end

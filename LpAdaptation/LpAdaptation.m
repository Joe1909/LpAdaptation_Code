function [out] = LpAdaptation(oracle, xstart, inopts)
%-------------------------------------------------------------------------
% Implementation of an algorithm for design centering and volume estimation, 
% proposal distribution is uniformly distributed
%
% Input:
% oracle: Name of the oracle function as string or function handle
% xstart: initial candidate solution, must be a feasible point
% inopts: option structure that determines internal strategy parameters (see code for details)
%
% Output:
% out: Output structure storing all relevant information (see code for
% details)
%
% When using this code please cite: 
% J. Asmus, C. L. Mueller and I. F. Sbalzarini. Lp-Adaptation: Simultaneous 
% Design Centering and Robustness Estimation of Electronic and Biological
% Systems. Scientific Reports 2017
%
% Josefine Asmus
% MOSAIC group, Center for Systems Biology Dresden (CSBD), Dresden, Germany
%
% Christian L. Mueller
% Flatiron Institute, Simons Foundation, New York City, NY, USA
%-------------------------------------------------------------------------

%% ---------- Returning default options when no argument are given ----------

if nargin==0
    
    N = 10
    disp('Default options for a N=10 dimensional problem')
    out = getDefaultOptions(N);
    return
    
end

%% ---------------------- Handling Input Parameters ----------------------

% ----- Handle oracle param
if isempty(oracle)
    error('Oracle not determined');
elseif ~ischar(oracle) && ~isa(oracle, 'function_handle')
    error('First argument ORACLE must be a string or function handle');
end
% If given oracle was a string convert it to handle and create universal handle for oracle
if ~isa(oracle, 'function_handle')
    oracle = str2func(oracle);
end

if (nargin==3) && isfield(inopts,'oracleInopts')
    oracleHdl = @(x) oracle(x, inopts.oracleInopts{:});
else
    oracleHdl = oracle;
end

% ----- Handle xstart param
if nargin < 2
    xstart = [];
end
if isempty(xstart)
    error('Initial search point, and problem dimension, not determined');
end
% dimension of the problem
N = length(xstart);

% ----- Handle inopts param. Merge options inopts and defopts
defopts = getDefaultOptions(N);
if nargin < 3 || isempty(inopts)
    inopts = [];
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

% Remember some other user's choices
isInitQprovided = isfield(inopts,'initQ');
isInitCprovided = isfield(inopts,'initC');
isbSavingOn = strcmp(opts.bSaving,'on');
isPlottingOn = strcmp(opts.Plotting,'on');

% Disable all 'lack of semicolon' warnings for clarity.
%#ok<*NOPRT>

%% ---------------------- Setup algorithmic Parameters ----------------------
% if it was specified if fixed schedule or no fixed schedule should be used
% hitting probability will be adapted (if not specified otherwise)
if isfield(inopts,'para_hitP_adapt') && isfield(inopts.para_hitP_adapt,'fixedSchedule') && ~isfield(inopts,'hitP_adapt')
    opts.hitP_adapt = 1;
end

maxMeanSize = opts.maxMeanSize;

if opts.hitP_adapt == 1 % adapting hitting probability
    valP = opts.para_hitP_adapt.PVec(1);
else
    valP = opts.valP;
end

% number of oracle output (first 0/1, then values that should be saved)
nOut = opts.nOut;

mueff=1; %#ok<NASGU> % ccov1 depends on mueff

% Algorithmic parameters
N_mu = myeval(opts.N_mu);   % weighting factor, controls how fast mean is shifted
N_C = opts.N_C;     % weights influence of accepted sample point on covariance matrix

ccov1 = myeval(opts.CMA.ccov1);

ccovmu      = myeval(opts.CMA.ccovmu);
PopSize     = myeval(opts.CMA.PopSize); %number of candidate solutions per iteration
cp          = opts.CMA.cp; %learning constant for evolution path
% expected length of input entries from sampling procedure
l_expected = sqrt(N);%get_expectedLength_Lp(N,pn,1e4); %dim, pnorm, numSamples

% windowSize for empirical hitting probability with moving window
windowSize = ceil(myeval(opts.windowSizeEval));%ceil(myeval(opts.windowSizeEval)/PopSize);
beta = myeval(opts.beta);
ss = myeval(opts.ss);       % expansion factor
sf = myeval(opts.sf);       % contraction factor
r = opts.r;         % radius (step size)
rMax = opts.MaxR;
rMin = opts.MinR;
condMax = opts.MaxCond; % max condition number
pn = opts.pn;       % norm of lp body
P_empAll = 0;
P_empWindow = 0;

if opts.hitP_adapt == 1 % if hitting probability should be adapted
    PVec        = opts.para_hitP_adapt.PVec;
    lPVec       = length(opts.para_hitP_adapt.PVec);
    
    if opts.para_hitP_adapt.fixedSchedule == 1 % if fixed schedule
        maxEvalSchedule = opts.para_hitP_adapt.maxEvalSchedule;
        numLastSchedule = opts.para_hitP_adapt.numLastSchedule;
        % check that PVec  has the same length as vectors maxEvalSchedule and numLastSchedule
        if length(PVec)~= length(maxEvalSchedule) || length(PVec)~= length(numLastSchedule) || length(maxEvalSchedule)~=length(numLastSchedule)
            error('vectors opts.para_hitP_adapt.PVec, opts.para_hitP_adapt.maxEvalSchedule and opts.para_hitP_adapt.numLastSchedule need to be of the same length');
        end
        %  values in maxEvalSchedule should sum up to 1
        if sum(maxEvalSchedule) < (1-1e-4) || sum(maxEvalSchedule) > (1+1e-4)
            warning('If values in maxEvalSchedule sum up to 1, you will test all hitting probabilites specified in opts.para_hitP_adapt.PVec and use allowed function evaluations')
            maxEvalSchedule
        end
        % all values in numLastSchedule should be in [0,1]
        if any(numLastSchedule<=0.1)||any(numLastSchedule >=0.9)
            warning('In numLastSchedule you can specify which portion of the samples should be used for calculation of final mean and radius. It should be in [0.1,0.9]') ;
            numLastSchedule(numLastSchedule<=0.1) = 0.1;
            numLastSchedule(numLastSchedule >= 0.9) = 0.9;
            numLastSchedule
        end
    else % no fixed schedule, adaption depending on step size, loewner and axis-alligned bound box volume approximation and hitting probability
        % check that meanOfLast is a number in [0,1]
        if opts.para_hitP_adapt.meanOfLast < 0 || opts.para_hitP_adapt.meanOfLast > 1
            warning(['opts.para_hitP_adapt.meanOfLast should be a number between 0 and 1. It is changed to ', num2str(defopts.para_hitP_adapt.meanOfLast)]);
            opts.para_hitP_adapt.meanOfLast = defopts.para_hitP_adapt.meanOfLast;
        end
        
        deviation_stepSize = opts.para_hitP_adapt.stepSize.deviation;
        deviation_VolApprox = opts.para_hitP_adapt.VolApprox.deviation;
        deviation_hitP = opts.para_hitP_adapt.hitP.deviation;
        deviation_stop = opts.para_hitP_adapt.deviation_stop;
        
        
        testEveryGen = ceil(myeval(opts.para_hitP_adapt.testEvery)/PopSize);
        testStartGen = ceil(myeval(opts.para_hitP_adapt.testStart)/PopSize);
        meanSize_stepSizeGen = ceil(myeval(opts.para_hitP_adapt.stepSize.meanSize)/PopSize);
        meanSize_VolApproxGen = ceil(myeval(opts.para_hitP_adapt.VolApprox.meanSize)/PopSize);
        meanSize_hitPGen = ceil(myeval(opts.para_hitP_adapt.hitP.meanSize)/PopSize);
        % check that testStartGen > 2* meanSizeGen
        if testStartGen < 2*max([meanSize_stepSizeGen,meanSize_VolApproxGen,meanSize_hitPGen])
            warning('opts.para_hitP_adapt.testStart needs to be at least 2 times bigger than opts.para_hitP_adapt....meanSize');
            testStartGen = 2* max([meanSize_stepSizeGen,meanSize_VolApproxGen,meanSize_hitPGen]);
            %error('opts.para_hitP_adapt.testStart needs to be at least 2 times bigger than opts.para_hitP_adapt....meanSize');
        end
        
    end % fixed schedule yes/no
else %hitP_adapt == 0
    numLast = myeval(opts.numLast);
    % check that numLast is not bigger than MaxEval
    if isfield(opts,'numLast') && numLast > opts.MaxEval
        numLast = myeval(defopts.numLast);
        warning(['Your chosen inopts.numLast is too big. It needs to be <= opts.MaxEval. It was changed to ',num2str(numLast)]);
    end
end %if hitting probability should be adapted

averageCovNum = myeval(opts.averageCovNum); % how many covariances used for averaged cov

% specified initial Covariance matrix or Cholesky matrix
% if nothing specified by user or if only Q is specified, use value in opts.initQ for Q
if ~isInitQprovided && ~isInitCprovided || isInitQprovided && ~isInitCprovided
    Q = opts.initQ;
    C = r^2*(Q*Q');
    % if Q and C are specified, check if that is consistent, use Q
elseif isInitQprovided && isInitCprovided
    Q = opts.initQ;
    if r^2*(Q*Q) == opts.initC
        C = opts.initC;
    else
        disp('specified Q and r yield this C: ');
        C = r^2*(Q*Q');
    end
    % if C is specified but not Q
elseif ~isInitQprovided && isInitCprovided
    C = opts.initC;
    % check if C is positiv semidefinit
    % symmetric (for real matrices) & all eigenvalues nonnegative
    [Bo,tmp]=eig(C);
    if (size(find(unique(C~=C')==1),1)>0) || any(diag(tmp<0))
        C
        diag(tmp)
        error('a covariance matrix needs to be positiv semidefinit!');
    elseif size(C,1)~=N %check if C and xstart have same size
        error('your input matrix C is not of the same size as your start vector');
    end
    
    diagD = diag(tmp);
    diagD = sqrt(diagD); % diagD contains standard deviations now
    detdiagD=(prod(diagD)); %normalize diagD
    diagD = diagD./(detdiagD.^(1/N));
    Q = Bo.*repmat(diagD',N,1); % O(n^2)
    r= nthroot(det(C),2*N);
end

%check if C and xstart have same size
if size(C,1)~=N
    error('your input matrix (Q or C) is not of the same size as your start vector');
end

[~,eigVals] = eig(Q*Q');
% Condition of initial C
condC = cond(Q*Q');

%% ----------------- Setup initial settings ---------------------------
% x_start needs to be a feasible point (inside body build up by the constraints)
% if point inside body -> c_T = 1, otherwise c_T = 0

% depending on the problem oracle tells if point is inside body or not
[outArgs{1:nOut}]= oracleHdl(xstart);
c_T = outArgs{1};

if c_T ~= 1
    error('x_start needs to be a feasible point!');
end

lastxAcc = xstart;

% Number of function evaluations
counteval = 1;
if opts.hitP_adapt == 1
    vcounteval=1;
 
    vcountgeneration=1;

    if opts.para_hitP_adapt.fixedSchedule == 1
        cntsave_Part=1;
    end
else
    % Number of evaluations after MaxEval - numLast evaluations
    countevalLast = 0;
    % Number of accepted points after MaxEval - numLast points
    lastNumAcc = 0;
end


countgeneration = 1;


% Number of all accepted points
% equals one because we need to start with a feasible point!
numAcc = 1;

% Number of accepted points for specific hitP
% equals one because we need to start with a feasible point!
vNumAcc = 1;
mu = xstart;

% ----------------- Setup output Parameters ---------------------------
if isbSavingOn
    % Trace of raw samples
    xRaw=nan(ceil(opts.MaxEval/opts.SavingModulo),N);
    xRaw(1,:)=xstart;
    
    % with accepted samples one can get upper bounds of the volume, also look
    % at samples mean
    xAcc=nan(opts.MaxEval,N); % all accepted x are saved
    cntAcc = nan(opts.MaxEval,1); % counteval of all accepted x
    cntAcc(1)=1;
    
    if nOut >1
        fxAcc = nan(opts.MaxEval,nOut-1);
        fxAcc(1,:) = outArgs{2:end};
    end
    xAcc(1,:)=xstart';
    
    if opts.unfeasibleSave ==1
        xNotAcc = nan(opts.MaxEval,N); % save all nonfeasible x as well
        cntNotAcc = nan(opts.MaxEval,1); %counteval of all infeasible x
        if nOut >1
            fxNotAcc = nan(opts.MaxEval,nOut-1);
        end
    end
    
    % Vector which tells if sample was accepted or not
    c_TVec=nan(ceil(opts.MaxEval/opts.SavingModulo),1);
    c_TVec(1)=c_T;
    
    if nOut >1
        fc_TVec = nan(ceil(opts.MaxEval/opts.SavingModulo),nOut-1);
        fc_TVec(1,:) = outArgs{2:end};
    end
    
    % Vector of evaluation indices when everything is saved
    countVec = nan(ceil(opts.MaxEval/opts.SavingModulo),1);
    countVec(1) = 1;
    
    % Cell array with stop flags
    stopFlagCell = cell(1,1);
    
    % Common settings for CMA/GaA
    VerboseModuloGen = ceil(opts.VerboseModulo/PopSize);
    SavingModuloGen  = ceil(opts.SavingModulo/PopSize);
    tmp_num = ceil(opts.MaxEval/SavingModuloGen);
    if opts.hitP_adapt == 1
        % count how often hitting probability is changed
        cntAdapt    = 1;
        % save different values for parameter when hitp changes
        N_muVec     = nan(lPVec,1);
        N_muVec(1)  = N_mu;
        betaVec     = nan(lPVec,1);
        betaVec(1)  = beta;            
        ssVec       = nan(lPVec,1);
        ssVec(1)    = ss;
        sfVec       = nan(lPVec,1);
        sfVec(1)    = sf;          
        iterVec     = nan(lPVec,1);
        windowSize = ceil(myeval(opts.windowSizeEval)/PopSize);
    end
    
    
    alpha_p = 1;
    pc = zeros(N,1); %initial evolution path for C
    %save generation of accepted samples
    cntAccGen = nan(opts.MaxEval,1);
    % Vector of evaluation indices when r,mu is saved
    cntVec = nan(tmp_num,1);
    cntVec(1) = 1;
    %do I need this?
    saveIndGeneration =2;

    if opts.hitP_adapt == 1            
        ccov1_Vec   = nan(lPVec,1);
        ccovmu_Vec  = nan(lPVec,1);
        cntGenVec   = nan(lPVec,1);
        PopSizeVec  = nan(lPVec,1);
        ccov1_Vec(1)   = ccov1;
        ccovmu_Vec(1)  = ccovmu;
        PopSizeVec(1) = PopSize;
    end
 
       
    % Vector of step lengths
    rVec=zeros(tmp_num,1);
    rVec(1)=r;
    
    % Vector of mu
    muVec = nan(tmp_num,N);
    muVec(1,:)=mu;
    
    % Vector of Volumina (of the lpBalls)
    volVec = zeros(tmp_num,1);
    volVec(1)=abs(det(Q))*Vol_lp(N,r,pn);
    
    % Vector of empirical acceptance probability
    P_empVecAll = zeros(tmp_num,1);
    P_empVecAll(1)=valP;
    P_empVecWindow = zeros(tmp_num,1);
    P_empVecWindow(1) = valP;
    numMuVec = zeros(windowSize,1);
    
    % Cell array of Q matrices
    if opts.bSaveCov==1
        QCell = cell(tmp_num,1);
        QCell{1} = Q;
    end
    
   
    
    if opts.hitP_adapt == 1
        % for each hitP save mu, r, hitP, Q
        muLastVec = nan(lPVec,N);
        rLastVec  = nan(lPVec,1);
        QLastCell = cell(lPVec,1);

        hitPLastVec = nan(lPVec,1);
        sizeLastVec = nan(lPVec,1);
        approxVolVec = nan(lPVec,1);
        
        if opts.para_hitP_adapt.fixedSchedule == 1
            MaxEval_Part = floor(maxEvalSchedule(cntAdapt)*opts.MaxEval);
            %when adapt == CMA, mu,r,Q,P needs only to be saved once for one generation
            numLast_Part = ceil(numLastSchedule(cntAdapt)*MaxEval_Part/PopSize);
            
            muLast = nan(numLast_Part,N);
            rLast  = nan(numLast_Part,1);
            P_empLast  = nan(numLast_Part,1);
   
        end
    else % if no adaptation of hitting probability
        %Trace of numLast mu, r and hitting probability values
        muLast      = nan(ceil(numLast/PopSize),N);
        rLast       = zeros(ceil(numLast/PopSize),1);
        P_empLast   = nan(ceil(numLast/PopSize),1);
        % Cell array of numLast Q matrices

    end
end %isbSavingOn

if opts.hitP_adapt == 1 %&& opts.para_hitP_adapt.fixedSchedule == 0
    % save all step sizes
    rVec_all = nan(ceil(opts.MaxEval/PopSize),1);
    rVec_all(1) = r;
    % save all hitP
    hitP_all = nan(ceil(opts.MaxEval/PopSize),1);
    hitP_all(1) = valP;
    % save corresponding function evaluations
    cnt_all= nan(ceil(opts.MaxEval/PopSize),1);
    cnt_all(1)=1;
end

% variable becomes 1 if adaptation of hitting probability occurs
van = 0;
% Final initialization
saveInd = 2;
saveIndAcc = 2;
saveIndNotAcc = 1; %everything starts with feasible point
stopFlag='';
if ~(opts.hitP_adapt == 1 && opts.para_hitP_adapt.fixedSchedule == 0)
    saveIndLast = 0;% Number of iterations after MaxEval - numLast evaluations
end
if opts.hitP_adapt == 1 && opts.para_hitP_adapt.fixedSchedule == 1
    MaxEval_Part = floor(maxEvalSchedule(cntAdapt)*opts.MaxEval);
    numLast_Part = floor(numLastSchedule(cntAdapt)*MaxEval_Part);
end

%% -------------------- Generation Loop --------------------------------

% several candidate solutions, adaptation of C as Adaptive Encoding, Hansen 2008
    PopSizeOld=[];
    lastEval = 1;

    [Bo,tmp]=eig(C);
    diagD=sqrt(diag(tmp));
    
    invB=diag(1./diagD)*Bo';
    while counteval(end) < (opts.MaxEval-PopSize-lastEval)%-lastEval
        counteval = counteval(end)+ (1:PopSize);
        if opts.hitP_adapt == 1
            vcounteval = vcounteval(end)+ (1:PopSize);
            vcountgeneration = vcountgeneration+1;
        end
        countgeneration = countgeneration+1;
        
        %% Generate PopSize candidate solutions uniformly distributed in lpBall
        % if norm >100 use infinity norm
        if pn>100 %sample uniform from hypercube with radius 1
            arz= -1 + 2.*rand(N,PopSize);%arz= -1 + 2.*rand(N,1);
        else
            %PopSize mutation vectors uniformly distributed in unit lp-ball in N dimensions
            %and with norm pn
            arz=lpBall(N,pn,PopSize)';
        end
        
        if any(isnan(Q(:)))
            disp('Q is nan');
            Q
        end
        
        arx = repmat(mu,1,PopSize) + r * (Q * arz);
        
        if isPlottingOn && isbSavingOn
            plotData(cntVec, countgeneration, saveIndGeneration, muVec, rVec, VerboseModuloGen, N, eigVals, r, P_empVecAll, P_empVecWindow, P_empAll, P_empWindow);
        end
        
        %% ORACLE
        % oracle gives 1 if point is in feasible region, 0 if not
        c_T = nan(PopSize,1);
        outArgsMat = nan(PopSize,nOut-1);
        
        for s=1:PopSize
            [outArgs{1:nOut}]= oracleHdl(arx(:,s));
            c_T(s) =outArgs{1};
            if nOut >1
                if any(cellfun('isempty',outArgs(2:end)))
                    outArgsMat(s,:) = NaN;
                else
                    outArgsMat(s,:) = outArgs{2:end};
                end
            end
        end
        
        % numfeas candidate solutions are in feasible region
        numfeas = sum(c_T==1);
        numMuVec(rem(countgeneration,windowSize)+1) = numfeas;
        numAccWindow = sum(numMuVec);
        
        if numfeas > 0
            % in pop all candidate solutions which are in feasible region
            pop = arx(:,c_T==1);
            weights = ones(numfeas,1)/numfeas; %uniform weights
            
            % count accepted solutions
            numAcc = numAcc + numfeas;
            if opts.hitP_adapt == 1
                vNumAcc = vNumAcc + numfeas;
            end
        end
        
        if opts.hitP_adapt == 1
            rVec_all(countgeneration)=r; %because counteval(1) is for initial values
            if van == 1
                P_empAll = valP;
                P_empWindow = valP;
                numMuVec = zeros(windowSize,1);
                numMuVec(rem(vcountgeneration,windowSize)+1) = numfeas;
                van = 0;
            else
                P_empAll = vNumAcc/vcounteval(end);
                cntEvalWindow = min(vcounteval(end),windowSize*PopSize);
                P_empWindow = numAccWindow/cntEvalWindow;
            end
            hitP_all(countgeneration) = P_empWindow;
            cnt_all(countgeneration)=counteval(end);

            if opts.para_hitP_adapt.fixedSchedule == 1
                % save muLast_Part mu and r values to get an average later
                if vcounteval(end) > (MaxEval_Part - numLast_Part)
                    muLast(cntsave_Part,:) = mu;
                    rLast(cntsave_Part) = r;
                    P_empLast(cntsave_Part) = P_empWindow;%P_empAll;
                  
                    
                    % if number of dedicated samples is spend change hitP
                    if vcounteval(end) >= (MaxEval_Part-PopSize+1)
                        disp('changing pVal');
                        muLastVec(cntAdapt,:) = mean(muLast(1:cntsave_Part,:));
                        rLastVec(cntAdapt)  = mean(rLast(1:cntsave_Part));
                        QLastCell{cntAdapt} = Q;
                        iterVec(cntAdapt)=counteval(end);
                        cntGenVec(cntAdapt) = countgeneration;
                        
                        hitPLastVec(cntAdapt) = P_empWindow;%P_empAll;
                        sizeLastVec(cntAdapt) = floor(numLast_Part/PopSize);
                        approxVolVec(cntAdapt) = P_empWindow * Vol_lp(N,rLastVec(cntAdapt),pn);%P_empAll * Vol_lp(N,rLastVec(cntAdapt),pn);
                        
                        if cntAdapt < lPVec
                            % restart with new hitP
                            cntAdapt = cntAdapt + 1;
                            
                            valP = PVec(cntAdapt);
                            PopSizeOld = PopSize;
                            PopSize=myeval(opts.CMA.PopSize);
                            if (mod(countgeneration,SavingModuloGen)~=0)
                                PopSizeOld = PopSize;
                            end
                            pc = zeros(N,1);%initial evolution path for C
                            mueff = 1;  %#ok<NASGU>
                            ccov1 = myeval(opts.CMA.ccov1);
                            ccovmu = myeval(opts.CMA.ccovmu);
                            beta = myeval(opts.beta);
                            ss = myeval(opts.ss);
                            sf = myeval(opts.sf);
                            N_mu = myeval(opts.N_mu);
                            
                           
                            
                            ccov1_Vec(cntAdapt) = ccov1;
                            ccovmu_Vec(cntAdapt) = ccovmu;
                            betaVec(cntAdapt)=beta;
                            ssVec(cntAdapt)=ss;
                            sfVec(cntAdapt)=sf;
                            N_muVec(cntAdapt)=N_mu;
                            PopSizeVec(cntAdapt) = PopSize;
                            windowSize = ceil(myeval(opts.windowSizeEval)/PopSize);
                            
                            MaxEval_Part = floor(maxEvalSchedule(cntAdapt)*opts.MaxEval);
                            numLast_Part = ceil(numLastSchedule(cntAdapt)*MaxEval_Part);
                            
                            muLast = nan(ceil(numLast_Part/PopSize),N);
                            rLast  = nan(ceil(numLast_Part/PopSize),1);
                            P_empLast  = nan(ceil(numLast_Part/PopSize),1);
                            cntsave_Part = 1;
                            
                            vcounteval = 1;
                            vcountgeneration = 1;
                            vNumAcc = 1;
                            van = 1;
                        else
                            stopFlag='all hitP (from PVec) tested';
                            break;
                        end %restart with new hitP
                    else % not all function evaluations spend yet (for this hitP)
                        cntsave_Part = cntsave_Part+1;
                    end
                end
            else %no fixed schedule for adaptation of hitting probability
                %save mu,r to get an average later
                muLast(vcountgeneration,:) = mu;
                rLast(vcountgeneration) = r;
                P_empLast(vcountgeneration) = P_empWindow;%P_empAll;
          
                % check after testStart iterations and then every testEvery iteration and save if last generation reached
                if ((vcountgeneration > testStartGen) && (mod(countgeneration,testEveryGen)==0)) || counteval(end) >= opts.MaxEval-PopSize-lastEval
                    mean1 = mean(rVec_all((countgeneration-2*meanSize_stepSizeGen+1):(countgeneration-meanSize_stepSizeGen)));
                    mean2 = mean(rVec_all((countgeneration-1*meanSize_stepSizeGen+1):(countgeneration)));
                    % also check hitting probability
                    mean3 = mean(hitP_all((countgeneration-2*meanSize_hitPGen+1):(countgeneration-meanSize_hitPGen)));
                    mean4 = mean(hitP_all((countgeneration-1*meanSize_hitPGen+1):(countgeneration)));
                    if (abs(mean1-mean2)/(mean1/2+mean2/2) < deviation_stepSize && ...
                            abs(mean3-mean4)/(mean3/2+mean4/2) < deviation_hitP) || counteval(end) >= opts.MaxEval-PopSize-lastEval
                        % check if volume approximation allows for
                        % adaptation of hitting probability
                        iidx1 = find(cntAccGen <= (countgeneration-meanSize_VolApproxGen),1,'last');
                        iidx2 = find(cntAccGen <= (countgeneration-1),1,'last');
                        if isempty(iidx1)
                            iidx1 = 1;
                        end
                        if isempty(iidx2)
                            iidx2 = 1;
                        end
                        x1 = xAcc(1:iidx1,:);
                        x2 = xAcc(1:iidx2,:);
                        
                        if size(x1,1) > 1 && size(x2,1) >1
                            epsLJ = 1e1;
                            %loewner of x1
                            [E,~] = lowner(x1',epsLJ);
                            C_loewner=inv(E);
                            Vol_loewner_x1 = Vol_lp(N,1,2) * sqrt(abs(det(C_loewner)));
                            %loewner of x2
                            [E,~] = lowner(x2',epsLJ);
                            C_loewner=inv(E);
                            Vol_loewner_x2 = Vol_lp(N,1,2) * sqrt(abs(det(C_loewner)));
                            
                            %axis alligned bounding box x1
                            mini = min(x1);
                            maxi = max(x1);
                            Vol_bb_x1 = prod(abs(maxi-mini));
                            %axis alligned bounding box x2
                            mini = min(x2);
                            maxi = max(x2);
                            Vol_bb_x2 = prod(abs(maxi-mini));
                            
                            if (abs(Vol_loewner_x1 - Vol_loewner_x2)/(Vol_loewner_x1/2 + Vol_loewner_x2/2) < deviation_VolApprox ...
                                    && abs(Vol_bb_x1-Vol_bb_x2)/(Vol_bb_x1/2+Vol_bb_x2/2) < deviation_VolApprox) || counteval(end) >= opts.MaxEval-PopSize-lastEval
                                % save r, Q, C, hitP
                                if cntAdapt >1
                                    numLastVec = ceil(cntGenVec(cntAdapt-1)+(vcountgeneration-cntGenVec(cntAdapt-1))*(1-opts.para_hitP_adapt.meanOfLast)):vcountgeneration;
                                else
                                    numLastVec = vcountgeneration - floor(vcountgeneration * opts.para_hitP_adapt.meanOfLast)+1:vcountgeneration;
                                end
                                muLastVec(cntAdapt,:) = mean(muLast(numLastVec,:));
                                rLastVec(cntAdapt)  = mean(rLast(numLastVec));
                                QLastCell{cntAdapt} = Q;
                                iterVec(cntAdapt)=counteval(end);
                                cntGenVec(cntAdapt) = countgeneration;
                                
                                hitPLastVec(cntAdapt) = P_empWindow;
                                sizeLastVec(cntAdapt) = vcountgeneration;
                                approxVolVec(cntAdapt) = P_empWindow * Vol_lp(N,rLastVec(cntAdapt),pn);%P_empAll * Vol_lp(N,rLastVec(cntAdapt),pn);
                               
                                
                                if cntAdapt < lPVec && counteval(end) < opts.MaxEval-PopSize-lastEval
                                    if cntAdapt >=2 && ...
                                            abs(approxVolVec(cntAdapt-1)-approxVolVec(cntAdapt))/(approxVolVec(cntAdapt-1)/2+approxVolVec(cntAdapt)/2) < deviation_stop
                                        disp('break loop because approximated volume does not change anymore');
                                        cntAdapt
                                        lPVec
                                        counteval
                                        opts.MaxEval
                                        lastEval
                                        stopFlag='no change in approxVol anymore';
                                        break;
                                    else
                                        disp('---------');
                                        disp('changing pVal');
                                        disp('---------');
                                        % restart with new hitP
                                        cntAdapt = cntAdapt + 1;
                                        valP = PVec(cntAdapt);
                                        PopSizeOld = PopSize;
                                        PopSize=myeval(opts.CMA.PopSize);
                                        if (mod(countgeneration,SavingModuloGen)~=0)
                                            PopSizeOld = PopSize;
                                        end
                                        meanSize_stepSizeGen = ceil(myeval(opts.para_hitP_adapt.stepSize.meanSize)/PopSize);
                                        meanSize_VolApproxGen = ceil(myeval(opts.para_hitP_adapt.VolApprox.meanSize)/PopSize);
                                        meanSize_hitPGen = ceil(myeval(opts.para_hitP_adapt.hitP.meanSize)/PopSize);
                                        testEveryGen = ceil(myeval(opts.para_hitP_adapt.testEvery)/PopSize);
                                        testStartGen = ceil(myeval(opts.para_hitP_adapt.testStart)/PopSize);
                                        
                                        if testStartGen < 2*max([meanSize_stepSizeGen,meanSize_VolApproxGen,meanSize_hitPGen])
                                            testStartGen = 2* max([meanSize_stepSizeGen,meanSize_VolApproxGen,meanSize_hitPGen]);
                                        end
                                        disp(testStartGen);
                                        % reinitialize
                                        mueff = 1; %#ok<NASGU>
                                        pc = zeros(N,1);%initial evolution path for C
                                        %evaluate parameters that depend on mueff and pVal
                                        ccov1 = myeval(opts.CMA.ccov1);
                                        ccovmu = myeval(opts.CMA.ccovmu);
                                        
                                        beta = myeval(opts.beta);
                                        ss = myeval(opts.ss);       % expansion factor
                                        sf = myeval(opts.sf);       % contraction factor
                                        N_mu = myeval(opts.N_mu);
                                        N_C = myeval(N_C);
                                        
                                        
                                        ccov1_Vec(cntAdapt) = ccov1;
                                        ccovmu_Vec(cntAdapt) = ccovmu;
                                        N_muVec(cntAdapt)=N_mu;
                                        betaVec(cntAdapt)=beta;
                                        ssVec(cntAdapt)=ss;
                                        sfVec(cntAdapt)=sf;
                                        PopSizeVec(cntAdapt) = PopSize;
                                        
                                        vcounteval = 1;
                                        vNumAcc = 1;
                                        van = 1;
                                        windowSize = ceil(myeval(opts.windowSizeEval)/PopSize);
                                    end
                                else
                                    disp('break loop because ~ cntAdapt < lPVec && counteval < opts.MaxEval-PopSize-lastEval');
                                    cntAdapt
                                    lPVec
                                    counteval
                                    opts.MaxEval
                                    lastEval
                                    PopSize
                                    stopFlag='maxcounteval or all hitP (from PVec) tested';
                                    break;
                                end
                            end %check volume approximation
                        end
                    end %check step size and hitP
                end % only check after testStart function evaluations and after every tetsEvery function evaluations
            end % end fixed schedule yes/no
        else % no adaptation of hitting probability
            cntEvalWindow = min(counteval(end),windowSize*PopSize);
            P_empAll = numAcc/counteval(end);
            P_empWindow = numAccWindow/cntEvalWindow;
            
            if (counteval(end)>=(opts.MaxEval-numLast))
                lastNumAcc = lastNumAcc+numfeas;
                saveIndLast = saveIndLast+1;
                countevalLast=countevalLast+PopSize;
                
                muLast(saveIndLast,:) = mu;
                rLast(saveIndLast) = r;
                P_empLast(saveIndLast) = P_empWindow;%P_empAll;
            
            end
        end % end adaptation of hitting probability yes/no
        
        if van == 0
            mueff = numfeas; %#ok<NASGU>
            %evaluate parameters that depend on mueff
            ccov1 = myeval(opts.CMA.ccov1);
            ccovmu = myeval(opts.CMA.ccovmu);
            N_mu = myeval(opts.N_mu);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Code for adaptation of C based on original Adaptive Encoding
        %    procedure by N. Hansen, see
        %    REFERENCE: Hansen, N. (2008). Adaptive Encoding: How to Render
        %    Search Coordinate System Invariant. In Rudolph et al. (eds.)
        %    Parallel Problem Solving from Nature, PPSN X,
        %    Proceedings, Springer. http://hal.inria.fr/inria-00287351/en/
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Adapt step size
        %%%%%%%%%%%%%%%%%%%%%%
        % adapt stepsize r depending on how many points were in feasible
        % region and how many were not, use parameters sf and ss as in GaA
        r = r*ss^numfeas*sf^(PopSize-numfeas);
        % Upper bound for step size to avoid numerical errors
        r = max(min(r,rMax),rMin);
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Adapt mu
        %%%%%%%%%%%%%%%%%%%%%%
        mu_old=mu;
        if numfeas ~= 0 && opts.madapt == 1
            mu = (1-1/N_mu)*mu + 1/N_mu*pop*weights;
            
            %%%%%%%%%%%%%%%%%%%%%
            % Update evolution paths
            %%%%%%%%%%%%%%%%%%%%%
            %adapt the encoding see Code Algorithm 7 in 'Adaptive Encoding for
            %Optimization' Nikolaus Hansen, 2008
            if (sum((invB*(mu - mu_old)).^2)==0)
                z=0;
            else
                alpha0 = l_expected/sqrt(sum((invB*(mu - mu_old)).^2));
                z= alpha0 *(mu - mu_old);
            end
        else
            z = 0;
        end
        
        pc = (1-cp)*pc+sqrt(cp*(2-cp))*z;
        S = pc*pc';
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Adapt covariance
        %%%%%%%%%%%%%%%%%%%%%%
        if numfeas ~=0 && opts.Cadapt == 1% no adapdation of C if no feasible solution in population
            
            arnorm = sqrt(sum((invB*(pop - repmat(mu_old,1,numfeas))).^2,1));
            alphai = l_expected * min(1./median(arnorm), 2./arnorm);
            
            if  opts.madapt == 0 || ~isfinite(alpha0)
                alpha0=1;
            end
            alphai(~isfinite(alphai)) = 1;
            
            zmu =  repmat(alphai, N,1) .* (pop - repmat(mu_old, 1, numfeas));
            Cmu = zmu * diag(weights) * zmu';
            
            C = (1-ccov1-ccovmu)*C+ccov1*alpha_p*S+ccovmu*(Cmu);
            C = triu(C)+triu(C,1)'; %enforce symmetry to prevent complex numbers
            
            [Bo, EV] = eig(C);
            EV = diag(EV);
            [EV,idx] = sort(EV); % vector of eigenvalues
            diagD = sqrt(EV);
            Bo = Bo(:,idx);
            
            if any(~isfinite(diagD))
                save(['tmp' opts.SaveFilename]);
                error(['function eig returned non-finited eigenvalues, cond(C)=' ...
                    num2str(cond(C)) ]);
            end
            if any(any(~isfinite(Bo)))
                error(['function eig returned non-finited eigenvectors, cond(C)=' ...
                    num2str(cond(C)) ]);
            end
            
            if ((condC < condMax) && condC >0)
                detdiagD=(prod(diagD)); %normalize D
                diagD = diagD./detdiagD.^(1/N);
                Q = Bo*diag(diagD);%faster
                invB = diag(1./diagD)*Bo';
            else
                if (mod(counteval,opts.VerboseModulo)==0)
                    disp( '-------------------------------------------');
                    disp(' Condition of C is too large and regularized');
                    disp( '-------------------------------------------');
                end
                % Regularize Q
                Q = Q + 1/N*eye(N,N);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%
        % Update condition of C
        QQ = (Q*Q');
        try
            [eigVecs,eigVals] = eig(QQ);
        catch
            Q
            counteval
        end
        
        condC = max(diag(eigVals))/min(diag(eigVals));
        
        if any(diag(eigVals)<0)
            disp('---------');
            disp('---------');
            saveInd
            size(xRaw)
            if ~isempty(xRaw)
                xRaw(max(1,saveInd-10):saveInd-1,:)
            end
            eigVecs
            eigVals
            condC
            det(Q)
            numfeas
            QQ
            warning('C contains negative eigenvalues')
            
            eigVals(eigVals<0)=1e-3;
            Q = eigVecs'*eigVals*eigVecs
            C = r^2*(Q*Q')
            r
            alphai
            alpha0
            mu_old
            mu
            disp('---------');
            disp('---------');
        end
        
        % Save all! accepted points
        if isbSavingOn
            if numfeas>0
                xAcc(saveIndAcc:(saveIndAcc+numfeas-1),:)=pop';
                if nOut>1
                    fxAcc(saveIndAcc:(saveIndAcc+numfeas-1),:) = outArgsMat(c_T==1,:);
                end
                cntAcc(saveIndAcc:saveIndAcc+numfeas-1,1)=counteval(c_T==1);
                cntAccGen(saveIndAcc:saveIndAcc+numfeas-1,1)=countgeneration;
                saveIndAcc = saveIndAcc+numfeas;
            end
            if opts.unfeasibleSave ==1
                if size(arx,2)~= PopSize
                    numunfeas = PopSizeVec(cntAdapt-1)-numfeas;
                else
                    numunfeas = PopSize-numfeas;
                end
                if numunfeas > 0
                    xNotAcc(saveIndNotAcc:(saveIndNotAcc+numunfeas-1),:)=arx(:,c_T==0)';
                    cntNotAcc(saveIndNotAcc:(saveIndNotAcc+numunfeas-1),:)=counteval(c_T==0);
                    if nOut>1
                        fxNotAcc(saveIndNotAcc:(saveIndNotAcc+numunfeas-1),:) = outArgsMat(c_T==0,:);
                    end
                end
                saveIndNotAcc = saveIndNotAcc+numunfeas;
            end
            
            
        end
        
        % Save r, mu (and possible Q) only (max) once in each iteration
        % since all candidate solutions are sampled from the same
        % distributon
        if (isbSavingOn && (mod(countgeneration,SavingModuloGen)==0)) || (counteval(end) > (opts.MaxEval-PopSize-lastEval)) 
            
            rVec(saveIndGeneration)=r;
            muVec(saveIndGeneration,:)=mu;
            volVec(saveIndGeneration) = abs(det(Q))*Vol_lp(N,r,pn);
            cntVec(saveIndGeneration) = counteval(end);
            
            %save evaluation number for which r,mu (Q) are stored
            P_empVecAll(saveIndGeneration)=P_empAll;
            P_empVecWindow(saveIndGeneration) = P_empWindow;
            % Cell array of Q matrices if covariances are stored
            if opts.bSaveCov==1
                QCell{saveIndGeneration}=Q;
            end
            
            
            
            saveIndGeneration = saveIndGeneration +1;
            
            %save rest
            if isempty(PopSizeOld)
                c_TVec(saveInd:saveInd+PopSize-1)=c_T;
                if nOut>1
                    fc_TVec(saveInd:saveInd+PopSize-1,:) = outArgsMat;
                end
                xRaw(saveInd:saveInd+PopSize-1,:)=arx';
                countVec(saveInd:saveInd+PopSize-1)=counteval;
                saveInd = saveInd+PopSize;
            else %if popsize changed with changed hitting probability
                c_TVec(saveInd:saveInd+PopSizeOld-1)=c_T;
                if nOut>1
                    fc_TVec(saveInd:saveInd+PopSize-1,:) = outArgsMat;
                end
                xRaw(saveInd:saveInd+PopSizeOld-1,:)=arx';
                countVec(saveInd:saveInd+PopSizeOld-1)=counteval;
                saveInd = saveInd+PopSizeOld;
                PopSizeOld = [];
            end
        end
        
        % command line message
        if (mod(countgeneration,VerboseModuloGen)==0)
            disp( '-------------------------------------------');
            disp([' Number of iterations: ',num2str(counteval(end))]);
            disp([' P_accAll:             ',num2str(P_empWindow)]);%P_empAll
            disp([' Search radius:        ',num2str(r)]);
        end
    end %while loop CMA


% Final sample is reserved for the final mu vector
if isbSavingOn
    
    out.xRaw = xRaw(1:saveInd-1,:);
    out.c_TVec = c_TVec(1:saveInd-1);
    if nOut >1
        out.fc_TVec = fc_TVec(1:saveInd-1,:);
    end
    out.countVec=countVec(1:saveInd-1);
    
    out.cntAcc=cntAcc(1:saveIndAcc-1,:);
    out.xAcc=xAcc(1:saveIndAcc-1,:);
    out.cntAccGen = cntAccGen(1:saveIndAcc-1,:);
    
    if nOut>1
        out.fxAcc = fxAcc(1:saveIndAcc-1,:);
    end
    
    out.countevalLast = counteval(end);
    
    if opts.unfeasibleSave ==1
        out.xNotAcc = xNotAcc(1:saveIndNotAcc-1,:);
        out.cntNotAcc = cntNotAcc(1:saveIndNotAcc-1,:);
        if nOut>1
            out.fxNotAcc = fxNotAcc(1:saveIndNotAcc-1,:);
        end
    end
    
 
    save_tmp = saveIndGeneration-1;
    out.cntVec = cntVec(1:save_tmp);
    out.rVec = rVec(1:save_tmp);
    out.muVec = muVec(1:save_tmp,:);
    out.P_empVecAll = P_empVecAll(1:save_tmp);
    out.P_empVecWindow = P_empVecWindow(1:save_tmp);
    out.volVec=volVec(1:save_tmp);
    
    if opts.bSaveCov==1
        out.QCell = QCell(1:save_tmp,:);
    end
    
    out.stopFlagCell = stopFlagCell;
    
    if opts.hitP_adapt == 0
        r=mean(rLast(1:saveIndLast));
        mu=mean(muLast(1:saveIndLast,:),1)';

    else %if hitP_adapt ==1
        
        out.adaptation.ccovmu = ccovmu_Vec(1:cntAdapt);
        out.adaptation.cntGenVec = cntGenVec(1:cntAdapt);
        out.adaptation.ccov1 = ccov1_Vec(1:cntAdapt);
        out.adaptation.PopSizeVec = PopSizeVec(1:cntAdapt);
      
        out.adaptation.N_muVec = N_muVec(1:cntAdapt);
        out.adaptation.betaVec = betaVec(1:cntAdapt);
        out.adaptation.ssVec = ssVec(1:cntAdapt);
        out.adaptation.sfVec = sfVec(1:cntAdapt);
        
        out.adaptation.iterVec = iterVec(1:cntAdapt);
        
        out.adaptation.muLastVec = muLastVec(1:cntAdapt,:);
        out.adaptation.rLastVec = rLastVec(1:cntAdapt);
        out.adaptation.QLastCell = QLastCell(1:cntAdapt);
        out.adaptation.hitPLastVec = hitPLastVec(1:cntAdapt);
        out.adaptation.sizeLastVec = sizeLastVec(1:cntAdapt);
        out.adaptation.approxVolVec = approxVolVec(1:cntAdapt);
        
        out.adaptation.rVec_all = rVec_all;
        out.adaptation.hitP_all = hitP_all;
        out.adaptation.cnt_all = cnt_all;
    end
    
    if opts.LastSaveAll == 1 && opts.hitP_adapt == 0
        out.outLast = outLast;
    end
end
r
mu
%% use last function evaluation to check if mean(muLast) is in feasible region!
%% if region is nonconvex or has holes it could happen that it is not feasible

% oracle gives 1 if point is in feasible region, 0 if not
[outArgs{1:nOut}]= oracleHdl(mu);
c_T = outArgs{1};

if c_T~=1
    warning('lastMu is not in feasible region!');
end

%% -------------------- Set Output --------------------------------
out.stopFlag = stopFlag;
if opts.hitP_adapt == 0
    out.P_empAll = numAcc/counteval(end);
    
    Vtest=Vol_lp(N,r,pn)*out.P_empAll;
    out.VolTest=Vtest;
    out.lastQ = Q;
    out.lastR = r;
    out.lastMu = mu;
    
    
else
    if isbSavingOn
        out.lastMu  = out.adaptation.muLastVec;
        out.lastR =  out.adaptation.rLastVec;
    end
end
out.P_emp = P_empAll;
out.P_empWindow = P_empWindow;
out.opts = opts;
out.opts.oracle = oracle; 

% -------------------- Ending Message ----------------------------
disp(['          Acceptance probability: ' num2str(P_empWindow)]);
disp(['          Stop flag: ' stopFlag]);

end

%-------------------------------------------------------------------------
%%%           get average covariance of numLast covariances            %%%
% ------------------------------------------------------------------------
function [out_averageCov,warningCell]= averageCov(outLast,N,pn,numLast)

p_test = pn;
dim = N;

Vec_similar_previous = nan(numLast-1,dim);

% save average eigenvalues and eigenvectors
eigVal_m_previous = nan(dim,1);
eigVec_m_previous = nan(dim);

eigValAll_previous = nan(numLast,dim);
warningCell=cell(numLast,1);
warningCnt = 1;

if numLast == 1
    %   use last covariance
    Q = outLast.QCell{end};
    r = outLast.rVec(end);
    C = r^2*(Q*Q');
    
    Vec_similar_previous=[];
    C_previous_sym=C;
    C_previous_asym=[];
    eigValAll_previous=[];
    
    warningCell{warningCnt} = 'only one covariance';
else
    for k =1:numLast
        %get MC Covariance in each step
        Q = outLast.QCell{end-numLast+k};
        r = outLast.rVec(end-numLast+k);
        C = r^2*(Q*Q');
        
        %eigen decomposition
        [eigVec,eigVals] = eig(C);
        
        if k==1
            %to save average eigenvalues and eigenvectors
            eigVal_m_previous=diag(eigVals);
            eigVec_m_previous=eigVec;
            
            %eigenvalues and eigenvectors from previous iteration
            eigVal_previous = diag(eigVals);
            eigVec_previous = eigVec;
            eigValAll_previous(k,:) = eigVal_previous;
        else
            check_allUsed_previous = zeros(dim,1);
            order_previous = zeros(dim,1);
            
            %signs
            si = ones(dim,1);
            
            for v1=1:dim
                
                %for each eigenvector, find corresponding eigenvector from
                %previous iteration and from .
                eigVec1 = eigVec(:,v1);
                
                %get vectors that are most alike from previous iteration
                [maxi_previous,idx_previous]=max(abs(eigVec1'*eigVec_previous));
                check_allUsed_previous(idx_previous) = check_allUsed_previous(idx_previous)+1;
                order_previous(v1)=idx_previous;
                
                %add eigenvalues
                diag_eigVals = diag(eigVals);
                eigVal_m_previous(idx_previous) = eigVal_m_previous(idx_previous)+diag_eigVals(v1);
                % check that vectors are pointing in same direction
                projection_previous = eigVec1'*eigVec_previous(:,idx_previous);
                if projection_previous < 0
                    eigVec1_previous = eigVec1*(-1);
                    si(v1)=-1;
                else
                    eigVec1_previous = eigVec1;
                end
                % add eigenvectors
                eigVec_m_previous(:,idx_previous)=eigVec_m_previous(:,idx_previous)+eigVec1_previous;
                
                %vec similar
                Vec_similar_previous(k-1,idx_previous) = maxi_previous;
            end
            if ~isempty(find(check_allUsed_previous==0, 1))% any(iszero(check_allUsed_previous))
                k
                check_allUsed_previous'
                error('not all eigenvectors matched to');
            end
            eigVal_previous = eigVals(order_previous);
            eigValAll_previous(k,:) = eigVal_previous;
            eigVec_previous = eigVec(:,order_previous);
            
            for v1=1:dim
                if si(v1)==-1
                    eigVec_previous(:,v1)= -1*eigVec_previous(:,v1);
                end
            end
        end
    end
    
    %divide by numLast --> get mean
    eigVec_previous = eigVec_m_previous./numLast;
    eigVals_previous = eigVal_m_previous./numLast;
    %normalize eigvectors v/norm(v)
    for d=1:dim
        eigVec_previous(:,d) = eigVec_previous(:,d)./norm(eigVec_previous(:,d)); % for k=1:dim eigVec2(:,k) = eigVec(:,k)./norm(eigVec(:,k));end
    end
    
    %get covariance
    C_previous = eigVec_previous *diag(eigVals_previous) * inv(eigVec_previous); %#ok<MINV>
    %symmetrize C
    C_previous_sym = 0.5*(C_previous + C_previous');
    % asymmetric part
    C_previous_asym = 0.5*(C_previous - C_previous');
end %if numLast == 1

out_averageCov.Vec_similar_previous = Vec_similar_previous;
out_averageCov.C_previous = C_previous_sym;
out_averageCov.C_previous_asym = C_previous_asym;
out_averageCov.p_test = p_test;
out_averageCov.dim = dim;
out_averageCov.eigValAll_previous = eigValAll_previous;

warningCell = warningCell(1:(warningCnt-1));
out_averageCov.warningCell = warningCell;
end

%-------------------------------------------------------------------------
%lower number of covariances used to average when error occurs in averageCov
function [C_p,C_normp] = run_averageCov(outLast,N,pn,numLast)
% if averageCov results in error, try again with smaller averageCovNum
count = 0;
err_count = 0;
while count == err_count
    try
        num = max(numLast - 2*err_count,1);% num = max(numLast-2^err_count+1,1);
        outCov = averageCov(outLast,N,pn,num);
        C_p = outCov.C_previous;
        % normalize C_p, such that det=1
        C_normp = C_p./det(C_p)^(1/N);
    catch
        err_count = err_count+1;
        num = max(numLast - 2*err_count,1); %num = max(numLast-2^err_count+1,1);
        disp(['try averageCov again with ',num2str(num), ' covariances']);
        if err_count > 10
            disp('try choosing inopts.averageCovNum smaller or MaxEval bigger such that the covariances that are used for averaging are more alike');
        end
    end
    count = count+1;
end
C_p = nearestSPD(C_p);
C_normp = nearestSPD(C_normp);
end

%-------------------------------------------------------------------------
%%%           From Mathworks, John D'Errico  30 Juli 2013         %%%
% ------------------------------------------------------------------------
function Ahat = nearestSPD(A)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

if nargin ~= 1
    error('Exactly one argument must be provided.')
end

% test for a square matrix A
[r,c] = size(A);
if r ~= c
    error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
    % A was scalar and non-positive, so just return eps
    Ahat = eps;
    return
end

% symmetrize A into B
B = (A + A')/2;

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[~,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
    [~,p] = chol(Ahat);
    k = k + 1;
    if p ~= 0
        % Ahat failed the chol test. It must have been just a hair off,
        % due to floating point trash, so it is simplest now just to
        % tweak by adding a tiny multiple of an identity matrix.
        mineig = min(eig(Ahat));
        Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
    end
end
end

% ---------------------------------------------------------------
% Volume of an Lp-ball
function [Vol] = Vol_lp(N,r,p)
% if p>100 we sample from a hypercube, so we also take the volume of a hypercube
% N - dimension
% r - radius
% p - p-norm
if p <= 100
    Vol = (2*gamma(1/p+1).*r).^N/(gamma(N/p+1));
else
    Vol = (2*r).^N;
end
end

% ---------------------------------------------------------------
function res=myeval(s)
if ischar(s)
    res = evalin('caller', s);
else
    res = s;
end
end

% ---------------------------------------------------------------
% Function to sample uniformly from Lp-ball
function[y]=lpBall(N,pn,number)
%-------------------------------------------------------------------------
% Generation of [number] samples from lpBall with norm [pn] in [n] dimension
% with center at [0 0 .. 0] and with radius 1
%
% Implementation of the Algorithm for real uniform generation
% see
% G. Calafiore, F. Dabbene, R. Tempo. Uniform Sample Generation in lp Balls for Probabilistic Robustness
% Analysis. Proceedings of the 37th IEEE Conference  on Decision & Control. Tampa, Florida USA. 1998
%
% Input:
% n:            Dimension
% pn:           p-norm (sum|xi|^p)^(1\p)
%   pn=1:       manhatten norm
%   pn=2:       euclidean norm
%   pn -> inf:  infinity norm
% number:       number of samples drawn from lp-ball
%
% Output:
% y:            real random vector y, which is uniformly distributed in
%               the lp-Ball B(1)
%-------------------------------------------------------------------------

%generate N independent random real scalars psi_i ~ G`(1/p,p) (generalized gamma distributed)
psi=gamrnd(1/pn,1,[number,N]);
psi=psi.^(1/pn);
%construct vector x e R^N of components xi=si*psi_i, where si are independent random signs
sign=randi(2,number,N); sign(sign==2)=-1;
X=psi.*sign;
%generate z=w^1/N, where w is random variable uniformly distributed in [0,1]
z=rand(number,1).^(1/N);
%y uniformly distributed in B(r)
y=repmat(z,1,N) .* ( X./repmat((sum(abs(X).^pn,2)).^(1/pn),1,N) );
end

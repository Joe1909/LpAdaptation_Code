clear inopts
clear out
addpath(genpath('../'));

oracle='oracleLpBall';
dim = 3;
Con = 3;

inopts.pn = 2; %p-norm of proposal distribution

%% set algorithm options
inopts.Plotting = 'on';%'off';

inopts.VerboseModulo=1e4;
inopts.MaxEval= dim*10^5;
inopts.SavingModulo =ceil(inopts.MaxEval/2000);

% adapt hitting probability with variable schedule
inopts.hitP_adapt = 1;
inopts.para_hitP_adapt.fixedSchedule = 0;
inopts.para_hitP_adapt.PVec = [0.35,0.15,0.06,0.03,0.01,0.005,0.002,0.001];
%inopts.para_hitP_adapt.maxEvalSchedule = [0.5,0.5];
inopts.para_hitP_adapt.testStart = 'max([2*opts.para_hitP_adapt.stepSize.meanSize,2*opts.para_hitP_adapt.hitP.meanSize,2*opts.para_hitP_adapt.VolApprox.meanSize])' ; 
% after how many iterations it should be tested for fist time; needs to be > 2*opts.para_hitP_adapt.meanSize
inopts.para_hitP_adapt.stepSize.meanSize = 'min(18/valP,maxMeanSize)';%'18/valP';%50; %compare two means with each other. meanSize defines the size of how how many r values are used for one mean
inopts.para_hitP_adapt.stepSize.deviation = 0.001; %if abs(mean1-mean2)/(mean1/2+mean2/2) < deviation, start step size adaptation
%loewner, axis-alligned bounding-box volume approximation
inopts.para_hitP_adapt.VolApprox.meanSize = 'min(18/valP,maxMeanSize)';
inopts.para_hitP_adapt.VolApprox.deviation = 0.001;
% hitting probability
inopts.para_hitP_adapt.hitP.meanSize = 'min(30/valP,maxMeanSize)';
inopts.para_hitP_adapt.hitP.deviation = 0.001;
% defining how many samples of each run are used for calculating the
% average if no fixed schedule for the changing of the hitting probability
inopts.para_hitP_adapt.meanOfLast = 1/4;%1/2; %number between 0 and 1; 
inopts.para_hitP_adapt.deviation_stop = 0.005;


%% feasible region
pn2 = 1; %p-norm
r2 = 1; %radius
mu2 = zeros(dim,1);
% square rooted covariance (C2 = Q2*Q2'*r2^2)
Q2  = diag(sqrt(logspace(0,Con,dim)));
Q2 = Q2./(det(Q2)^(1/dim));

inopts.oracleInopts{1} = pn2; 
inopts.oracleInopts{2} = r2; 
inopts.oracleInopts{3} = mu2;
inopts.oracleInopts{4} = Q2; 
 
% true volume
Vol_t = Vol_lp(dim,r2,pn2);
numRep = 2;                 
outCell = cell(numRep,1);

for rep=1:numRep
    disp(rep);
    tmp = lpBall(dim,pn2,1)';
    xstart = mu2 + r2*(Q2*tmp);
    
    out = LpAdaptation(oracle,xstart,inopts);
    outCell{rep} = out; 
end
volVec = nan(numRep,length(inopts.para_hitP_adapt.PVec));



figure; 
cmap=colormap(lines); 
for rep=1:numRep
    plot(outCell{rep}.cntVec, outCell{rep}.volVec.*outCell{rep}.P_empVecWindow,'col',cmap(rep,:));
    hold on
    %volVec(rep,:) = outCell{rep}.adaptation.approxVolVec; 
    plot([0 inopts.MaxEval],[Vol_t Vol_t],'k--','LineWidth',2);
    xlabel('evaluations');
    ylabel('estimated volume');
end
                    
                    
                    


                
     



clear inopts
clear out
addpath(genpath('../'));

oracle='oracle_packman';
xstart=[-1;0];

% dimension of the problem
dim = length(xstart);

%set up designCentering options
% p-norm of proposal distribution
inopts.pn = 2;

% number of maximal evaluations
inopts.MaxEval= 1e4;
% saving after every i-th generation
inopts.SavingModulo =10;
% how often output in command window and plot (if plot is on)
inopts.VerboseModulo=1000;
inopts.Plotting = 'on';

% adapt hitting probability
inopts.hitP_adapt = 1;

% how to adapt hitting probability
inopts.para_hitP_adapt.PVec = [0.35,0.55,0.75,0.9];
% use fixed schedule
inopts.para_hitP_adapt.fixedSchedule = 1; % 1: fixed schedule, 0: no fixed schedule, then it will be decided on the fly when to change the hitting probability
% of interest if hitP_adapt ==1 and fixed schedule
inopts.para_hitP_adapt.maxEvalSchedule = (1/4)*ones(1,4);%[1/2,1/8,1/8,1/8,1/8];%[3/5,3/50,4/50,6/50,7/50]; % proportion of maxEval 
% that are used for every run with a different hitting probability, 
% vector should have same length as PVec, sum of this vector should ideally
% equal 1 (otherwise function evaluations are wasted)
inopts.para_hitP_adapt.numLastSchedule = [0.5,0.5,0.5,0.5];%[1/2,3/4,3/4,3/4,3/4]; % over how many 
% samples of each run should be averaged to get radius r and mean mu
% of found feasible region


% start algorithm
out = LpAdaptation(oracle,xstart,inopts);

%% Figure
 cols =[   0.8941    0.1020    0.1098
    0.2157    0.4941    0.7216
    0.3020    0.6863    0.2902
    0.5961    0.3059    0.6392
    1.0000    0.4980         0
    1.0000    1.0000    0.2000
    0.6510    0.3373    0.1569
    0.9686    0.5059    0.7490
    0.6000    0.6000    0.6000
    0.3216    0.3216    0.3216
         0         0         0];

[X1,Y1] = meshgrid(linspace(-1,1,500),linspace(-1,1,500));
[a,b]=size(X1);
sample = [X1(:),Y1(:)];
inside1 = oracle_packman(sample);
Z1 = reshape(inside1,[a,b]);

iterSave = length(out.rVec);
pn = out.opts.pn;
cnt_old = 1;
 
f2=figure(2);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
[C H]=contourf(X1,Y1,Z1,[0 1]);
set(H,'LineWidth',2);
colormap([1 1 1;0.7 0.7 0.7]);
xlabel('x','FontSize',30)
ylabel('y','FontSize',30)
axis square%equal
grid on
xlim([-2,2]);
ylim([-2,2]);
    
    
for iter = [1:2:iterSave]

    hold on;
    eval_i = out.cntVec(iter);
    cnt = find(out.cntAcc <= eval_i,1,'last');
    plot(out.xAcc(cnt_old:cnt,1),out.xAcc(cnt_old:cnt,2),'r.');
    grid on
    cnt_old = cnt; 
    mu_tmp = out.muVec(iter,:);
    Q_tmp = out.QCell{iter};
    r_tmp = out.rVec(iter);
            
    ball=lpBall2surface(4000,2,pn,r_tmp,mu_tmp',Q_tmp,0);
    set(gca,'fontsize',18);
    
    idx_col = find(out.cntVec(iter)< out.adaptation.iterVec,1,'first')+1;
    
    h = plot(ball(:,1),ball(:,2),'.','col',cols(idx_col,:),'MarkerSize',6);
    m = plot(mu_tmp(1),mu_tmp(2),'x','col',cols(idx_col,:),'MarkerSize',15,'linew',3);
    set(gcf, 'Color', 'w');
    txt1 = ['P = ',num2str(out.opts.para_hitP_adapt.PVec(idx_col-1))];
    t = text(-3.7,3.7,txt1,'fontSize',30,'col',cols(idx_col,:));
    pause(0.01)
    
    if iter<iterSave
        delete(h);
        delete(m);
        delete(t);
    end
end
hold off





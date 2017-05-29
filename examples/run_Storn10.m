clear all 
close all
addpath(genpath('../'));

%volume of problem storn 1999 fig 10
V2=101;
%constraints see Storn Paper 1999 Fig.10
xConstraints = [0,5,5,10,10,5,5,0,0];
yConstraints = [-0.1,-0.1,-10,-10,10,10,0.1,0.1,-0.1];    

oracle='oracleStorn10';
xstart=[0;0];

% dimension of the problem
dim = length(xstart);
%set up designCentering options
inopts.pn = 2; %change p-norm
inopts.MaxEval= 2000;
inopts.SavingModulo =100;
inopts.VerboseModulo=100;

inopts.Plotting = 'on';

out = LpAdaptation(oracle,xstart,inopts);


f2=figure(2);

iterSave = length(out.rVec);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
plot(xConstraints,yConstraints);
[~,~,~,current_entries]=legend;
 legend([current_entries {'boundary of feasible region'}],'Location','NorthWest','fontsize',20)
xlim([-1 12]);
ylim([-15 15]);
hold on;
pn = out.opts.pn;
plot(out.xAcc(1,1),out.xAcc(1,2),'b.');
cnt_old = 1;
t2 = text(8,14,['evaluations 1'],'HorizontalAlignment','left','VerticalAlignment','top','fontsize',18);

for iter=1:1:iterSave
    idx = out.cntVec(iter);
    cnt = find(out.cntAcc <= idx,1,'last');
    plot(out.xAcc(cnt_old+1:cnt,1),out.xAcc(cnt_old+1:cnt,2),'b.');
    cnt_old = cnt; 
    mu = out.muVec(iter,:);
    Q = out.QCell{iter};
    r = out.rVec(iter);
    
    [eigvec,eigval] = eig(r^2*Q*Q');
    
    ball=lpBall2surface(300000,2,pn,r,mu',Q,0);
    
    h = plot(ball(:,1),ball(:,2),'r.','markerSize',3);
    m = plot(mu(1),mu(2),'r.');
    
    delete(t2);
    t2 = text(8,14,['evaluations ',num2str(idx)],'HorizontalAlignment','left','VerticalAlignment','top','fontsize',18);
    pause(0.1)
    if iter<iterSave
        delete(h);
        delete(m);
    end
    
    

end
    
    

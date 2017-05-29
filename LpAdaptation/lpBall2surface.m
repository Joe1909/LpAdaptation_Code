
function[y]=lpBall2surface(numSamples,N,pn,r,center,Q,t)
%-------------------------------------------------------------------------
% Implementation of the Algorithm for real uniform generation
% see 
% G. Calafiore, F. Dabbene, R. Tempo. Uniform Sample Generation in lp Balls for Probabilistic Robustness
% Analysis. Proceedings of the 37th IEEE Conference  on Decision & Control. Tampa, Florida USA. 1998
%
% Input:
% numSamples:   number of desired samples
% n:            Dimension
% pn:           p-norm (sum|xi|^p)^(1\p) 
%   pn=1:       manhatten norm
%   pn=2:       euclidean norm
%   pn -> inf:  infinity norm
% r:            radius of the ball          
% center:       center of the ball          
% Q:            square root of covariance matrix C=r^2Q*Q'
% t:            plot lpBall: t=1, otherwise t=0
%
% Output:
% y:            real random vectors y, which are uniformly distributed on
% the surface of the lp-Ball B(r)
%-------------------------------------------------------------------------

%generate n independent random real scalars psi_i ~ G`(1/p,p) (generalized
%gamma distributed)
psi=gamrnd(1/pn,1,[numSamples,N]);
psi=psi.^(1/pn);
%construct vector x e R^n of components xi=si*psi_i, where si are
%independent random signs
sign=randi(2,numSamples,N); sign(sign==2)=-1; 
X=psi.*sign;
%generate z=w^1/n, where w is random variable uniformly distributed in
%[0,1]
% % z=(rand(numSamples,1)).^(1/N);%z=(my_rand([numSamples,1])).^(1/N); %
% % %y uniformly distributed in B(r)
% % y=repmat(r*z,1,N).*(X./repmat((sum(abs(X).^pn,2)).^(1/pn),1,N));

%y uniformly distributed on border of B(r)
y=r*(X./repmat((sum(abs(X).^pn,2)).^(1/pn),1,N));
y=(Q*y')';

%move ball to center
y=y+repmat(center',numSamples,1);

%plot if t==1
if(t==1)
    if N>=3
        plot3(y(:,1),y(:,2),y(:,3),'r.','markersize',1); 
        axis equal;rotate3d off; rotate3d on;drawnow;shg; 
        xlabel('x1');ylabel('x2');zlabel('x3');
        title(sprintf('numSamples %d pnorm %4.2f',numSamples,pn));
    else
        plot(y(:,1),y(:,2),'r.','markersize',1); 
        %axis equal;zoom off; zoom on;drawnow;shg; 
    end
end
end

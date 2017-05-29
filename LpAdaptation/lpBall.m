
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
% pn:           p-norm (sum|xi|^p)^(1\p)    default: 2
%   pn=1:       manhatten norm
%   pn=2:       euclidean norm
%   pn -> inf:  infinity norm
% number:       number of samples drawn from lp-ball
%
% Output:
% y:            real random vector y, which is uniformly distributed in
%               the lp-Ball B(1)
%-------------------------------------------------------------------------

% % % % ---------------------- Handling Input Parameters ----------------------

if isempty(N)
    error('dimension not determined');
end

if nargin < 2
    pn = 2;
    number = 1;
end
if nargin < 3
    number = 1;
end

%generate N independent random real scalars psi_i ~ G`(1/p,p) (generalized
%gamma distributed)
psi=gamrnd(1/pn,1,[number,N]);
psi=psi.^(1/pn);
%construct vector x e R^N of components xi=si*psi_i, where si are
%independent random signs
sign=randi(2,number,N); sign(sign==2)=-1; 
X=psi.*sign;
%generate z=w^1/N, where w is random variable uniformly distributed in
%[0,1]
z=rand(number,1).^(1/N);%z=(my_rand([1,1])).^(1/n); %
%y uniformly distributed in B(r)
y=repmat(z,1,N).*(X./repmat((sum(abs(X).^pn,2)).^(1/pn),1,N));
end

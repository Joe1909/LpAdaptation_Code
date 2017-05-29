function [f]=oracleLpBall(x,np,r,mu,Q)
%%%% INPUT: %%%%
% x candidate solution(s)
%   single candidate solution as column vector
%   multiple candidate solutions: each column one candidate solution
% np p-norm of LpBall
% r radius of LpBall
% mu center of LpBall
% Q deformation matrix
%%%% OUTPUT: %%%%
% 1 if x inside LpBall
% 0 else
% when x was an array, the output is a vector of 0s and 1s

%%
if size(x,2)==1
    xtest = inv(Q)/r*(x-mu);
    
    if np > 100
        if max(abs(xtest))<=1
            f=1;
        else
            f=0;
        end
    else
        if ((sum(abs(xtest).^np))^(1/np))<=1
            f=1;
        else
            f=0;
        end
    end
else
    number = size(x,2);
    %     x
    %     mu
    %     number
    b=(x-repmat(mu,1,number));
    xtest = inv(Q)/r*(b);
    
    if np > 100
        f = max(abs(xtest))<=1;
    else
        f=(((sum(abs(xtest).^np)).^(1/np))<=1)';
    end
    
end
end



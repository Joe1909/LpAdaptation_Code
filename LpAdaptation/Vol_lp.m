function [Vol] = Vol_lp(N,r,p)
    % N - dimension
    % r - radius
    % p - p-norm 
    if p > 100
        Vol = (2*r).^N;
    else
        Vol = (2*gamma(1/p+1).*r).^N/(gamma(N/p+1));
    end
function inside = oracle_packman(x)

if isvector(x)
    n = length(x);
    if iscolumn(x)
        x = x';
    end
    number=1;
else
    [number,n]=size(x);
end


if n~=2
    error('Wrong dismension of input. Test body is in 2D')
end

mu = [0.3,0];
Q = [2.5,0;0,1];
r = 0.3; 

if number==1
    f1 = (((sum(abs(x).^2))^(1/2))<=1);
    f2 = (((sum(abs(inv(Q)/r*(x-mu)').^2))^(1/2))<=1);
    
else
    x1=x';
    f1=(((sum(abs(x1).^2)).^(1/2))<=1)';
    %     x
    %     mu
    %     number
    mu1=mu';
    b=(x1-repmat(mu1,1,number));
    xtest = inv(Q)/r*(b);
    
   
    f2=(((sum(abs(xtest).^2)).^(1/2))<=1)';
    
    
end

inside = f1 & ~f2;
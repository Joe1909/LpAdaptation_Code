function inside = oracle_heavyTailedStar(x)

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
a = 1.1;%1.5; 
b = 1.35; 
mu1 = [b,0];
Q1 = [a,0;0,1];
r1 = 0.9; 

mu2 = [-b,0];
Q2 = [a,0;0,1];

mu3 = [0,b];
Q3 = [1,0;0,a];

mu4 = [0,-b];
Q4 = [1,0;0,a];

if number==1
    f0 = (((sum(abs(0.5*x).^2))^(1/2))<=1);
    f1 = (((sum(abs(inv(Q1)/r1*(x-mu1)').^2))^(1/2))<=1);
    f2 = (((sum(abs(inv(Q2)/r1*(x-mu2)').^2))^(1/2))<=1);
    f3 = (((sum(abs(inv(Q3)/r1*(x-mu3)').^2))^(1/2))<=1);
    f4 = (((sum(abs(inv(Q4)/r1*(x-mu4)').^2))^(1/2))<=1);
else
    x1=x';
    f0=(((sum(abs(0.5*x1).^2)).^(1/2))<=1)';
    %     x
    %     mu
    %     number
    mu=mu1';
    b=(x1-repmat(mu,1,number));
    xtest = inv(Q1)/r1*(b);
    f1=(((sum(abs(xtest).^2)).^(1/2))<=1)';
    
    mu=mu2';
    b=(x1-repmat(mu,1,number));
    xtest = inv(Q2)/r1*(b);
    f2=(((sum(abs(xtest).^2)).^(1/2))<=1)';
    
    mu=mu3';
    b=(x1-repmat(mu,1,number));
    xtest = inv(Q3)/r1*(b);
    f3=(((sum(abs(xtest).^2)).^(1/2))<=1)';
    
    mu=mu4';
    b=(x1-repmat(mu,1,number));
    xtest = inv(Q4)/r1*(b);
    f4=(((sum(abs(xtest).^2)).^(1/2))<=1)';
    
    
end

inside = f0 & ~f1 &~f2 &~f3 &~f4;
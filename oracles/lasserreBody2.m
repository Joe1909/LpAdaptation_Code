function inside = lasserreBody2(x)

if isvector(x)
    n = length(x);
    if iscolumn(x)
        x = x';
    end
else
    [~,n]=size(x);
end

if n~=2
    error('Wrong dismension of input. TEs body is in 2D')
end

g = sum(x.^6,2)-1.925*prod(x.^3,2);
inside = (g<1);


function [f]=oracleStorn10(x)

if size(x,2)==1
    if x(1)>=0&&x(1)<=5&&x(2)<=0.1&&x(2)>=-0.1
        f=1;
    elseif x(1)>=5&&x(1)<=10&&x(2)<=10&&x(2)>=-10
        f=1;
    else
        f=0;
    end
else
   len = size(x,2);
    f=zeros(len,1);
    for k=1:len
        if x(1,k)>=0&&x(1,k)<=5&&x(2,k)<=0.1&&x(2,k)>=-0.1
            f(k)=1;
        elseif x(1,k)>=5&&x(1,k)<=10&&x(2,k)<=10&&x(2,k)>=-10
            f(k)=1;
        else
            f(k)=0;
        end
    end 
end
end
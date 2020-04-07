function [y2,x,r2] = truncNorm(par,x) 

r2 = round(max(0,par(1) + par(2) * randn(1e6,1)));

if nargin==1
    [y2,x]=hist(log10(r2),100);
else
    y2=hist(log10(r2),x);
end

end


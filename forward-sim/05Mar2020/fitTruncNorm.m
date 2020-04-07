
function L = fitTruncNorm(par,x, r1)
[y2,~,r2]=truncNorm(par,x);

L = 10*(mean(r1)-mean(r2)).^2+sum((cumsum(hist(r1,x))-cumsum(y2)).^2/sum(y2)^2);
end

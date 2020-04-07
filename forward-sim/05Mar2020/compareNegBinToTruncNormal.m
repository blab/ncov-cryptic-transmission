% negbinom vs truncated normal

clc; 
figure(4); clf;

r=0.54;
mu=3.2;

p = 1/(1+mu/r);

r1 = nbinrnd(r,p,1e6,1);
mean(r1)
[y,x]=hist(log10(r1),100);

clf; hold on;
plot(x,cumsum(y),'r')


m=-80*7/365;
s=400*7/365;

par0(1)=-700*7/365;
par0(2)=900*7/365;

fitPar=fminsearch( @(par) fitTruncNorm(par,x, r1),par0);

fitTruncNorm(fitPar,x, r1)
[y2,x,r2]=truncNorm(fitPar,x);

mean(r2)
plot(x,cumsum(y2),'b')
legend('negative binomial','truncated normal')

round(fitPar* 365/8)

% normal has a less fat right tail, but the concept is very similar.
% Honestly, I believe the trunc normal for human to human transmission, vs
% aerosol environmental.
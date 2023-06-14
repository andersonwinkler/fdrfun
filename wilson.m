function LU = wilson(n,X,alpha)
% Compute the confidence interval for a set of Bernoulli
% trials using the Wilson method.
%
% Usage:
% LU = confint(n,X,alpha)
%
% - n     : Number of trials.
% - X     : Number of successful trials.
% - alpha : Coverage of the confidence interval.
%           For 95% CI, use alpha = 0.05.
% - LU    : Two-element vector with lower and upper CI.
%
% The variables n, X and alpha can be either scalars
% or arrays. When using arrays, use consistent sizes.
%
% References:
% * Wilson EB. Probable inference, the law of succession, and
%   statistical inference. JASA. 1927 22(158):209-212.
% * Brown LD, Cai TT, DasGupta AA. Interval estimation for a
%   binomial proportion. Statistical Science. 2001 16(2):101-133.
%
% _____________________________________
% Anderson Winkler & Tom Nichols
% FMRIB / University of Oxford
% Apr/2012
% http://brainder.org

k  = norminv(1-alpha/2);
p  = X./n;
q  = 1 - p;
Xt = X + (k.^2)/2;
nt = n + k.^2;
pt = Xt./nt;
L = pt - k.*sqrt(n.*p.*q + (k.^2)/4)./nt;
U = pt + k.*sqrt(n.*p.*q + (k.^2)/4)./nt;
LU = [L U];

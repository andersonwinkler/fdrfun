function varargout = simes(varargin)
% Computes the Simes-adjusted p-value for a vector of p-values.
%
% Usage:
% [pthr,padj] = simes(pvals)
%               simes(pval,q)
%
% Inputs:
% pvals  = Vector of p-values.
% q      = Allowed proportion of false positives (q-value).
%          Default = 0.05.
%
% Outputs:
% pthr   = Simes threshold; if any pval is below, reject the global.
% padj   = Simes-adjusted p-values.; if any padj<=q, reject the global.
%
% Reference:
% * Simes RJ. An improved Bonferroni procedure for multiple tests
%   of significance. Biometrika. Biometrika Trust; 1986;73(3):751. 
%
% ________________________________
% Anderson M. Winkler
% Univ of Texas Rio Grande Valley
% Sep/2025
% http://brainder.org

% Accept arguments
switch nargin
    case 0
        error('Error: Not enough arguments.');
    case 1
        pval = varargin{1};
        qval = 0.05;
    case 2
        pval = varargin{1};
        qval = varargin{2};
    otherwise
        error('Error: Too many arguments.')
end

% Check if pval is a vector
if numel(pval) ~= length(pval)
    error('p-values should be a row or column vector, not an array.')
end

% Check if pvals are within the interval
if numel(pval)>0 && (min(pval) < 0 || max(pval) > 1)
    error('Values out of range (0-1).')
end

% Check if qval is within the interval
if qval < 0 || qval > 1
    error('q-value out of range (0-1).')
end

% ========[PART 1: SIMES THRESHOLD]========================================

% Sort p-values
[pval,~] = sort(pval);

% Number of observations
V = numel(pval);

% Order (indices), in the same size as the pvalues
idx = reshape(1:V,size(pval));

% Line to be used as cutoff
thrline = idx*qval/V;

% Find the largest pval, still under the line
thr = max(pval(pval<=thrline));

% Deal with the case when all the points under the line
% are equal to zero, and other points are above the line
if thr == 0
    thr = max(thrline(pval<=thrline));
end

% Case when it does not cross
if isempty(thr)
    thr = 0;
end

% Returns the result
varargout{1} = thr;

% ========[PART 2: SIMES ADJUSTED]========================================

if nargout == 2
    
    % p-adjusted
    padj = min(pval.*V./idx);
    varargout{2} = padj;
end

% That's it!

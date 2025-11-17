function [pthr,padj,qfac] = bb2014(varargin)
% Computes the FDR-threshold for multiple families of hypotheses
% using the Benjamini-Bogomolov selective inference procedure.
%
% Usage:
% [pthr,padj,qfac] = bb2014(Pset)
%                    bb2014(Pset,q)
%                    bb2014(Pset,q,fdrfun)
%                    bb2014(Pset,q,fdrfun,cV)
%
% Inputs:
% Pset   = Set of p-values organized as families. Can be:
%          - Cell array: each cell contains a vector of p-values for one family.
%          - Matrix: each row contains p-values for one family.
% q      = Allowed proportion of false positives (q-value).
%          Default = 0.05.
% fdrfun = Function handle for FDR method to use.
%          Default = @bky7.
% cV     = If set to anything but 1, uses an harmonic sum for c(V).
%          Default = 1.
%
% Outputs:
% pthr   = FDR thresholds for each family (cell array or matrix).
% padj   = FDR adjusted p-values for each family (cell array or matrix).
% qfac   = Ratio R/S, i.e., number of familes with significant results
%          after Simes divided by the number of families.
%          To threshold padj, use q*qfac.
% 
% The procedure works in three steps:
% 1. Applies Simes test to each family of hypotheses.
% 2. Applies FDR correction to the Simes p-values across families.
% 3. Applies FDR to each selected family with adjusted q-value (q*R/nFam).
%
% Families not selected in step 2 receive pthr = 0 and padj = 1.
% 
% Note that the adjusted p-values do **not** depend on the supplied q-value,
% but they do depend on the choice of c(V).
%
% References:
% * Benjamini Y, Bogomolov M. Selective inference on multiple
%   families of hypotheses. Journal of the Royal Statistical
%   Society: Series B (Statistical Methodology).
%   2014 Jan 18;76(1):297–318.
%
% ________________________________
% Anderson M. Winkler
% UTRGV
% Sep/2025
% http://brainder.org

narginchk(1,3);
Pset   = varargin{1};
qval   = 0.05;
fdrfun = @bky7;
cV     = 1;
if nargin > 1 && ~isempty(varargin{2})
    qval   = varargin{2};
end
if nargin > 2 && ~isempty(varargin{3})
    fdrfun = varargin{3};
end
if nargin > 3
    if varargin{4}
        cV = 1;
    else
        cV = sum(1./(1:numel(pval)));
    end
end

% ========[PART 1: SIMES ON EACH FAMILY]===================================
if iscell(Pset)
    nFam = length(Pset);
    if nFam ~= numel(Pset)
        error('Input cell arrays must have only 1 dimension')
    end
    psimes = zeros(nFam,1);
    for fam = 1:nFam
        if isempty(Pset{fam})
            psimes(fam) = 1;
        else
            [~,psimes(fam)] = simes(Pset{fam},qval);
        end
    end
elseif ismatrix(Pset)
    nFam = size(Pset,1);
    psimes = zeros(nFam,1);
    for fam = 1:nFam
        [~,psimes(fam)] = simes(Pset(fam,:),qval);
    end
else
    error('Input set of p-values must be a matrix or cell array');
end

% ========[PART 2: FDR ON SIMES]===========================================
[~,psimesfdradj] = fdrfun(psimes,qval,cV); % p after FDR on Simes
Ridx = psimesfdradj <= qval;
R    = sum(Ridx);
qfac = R/nFam;

% ========[PART 3: FDR ON EACH SELECTED FAMILY]============================
if iscell(Pset)
    pthr = cell(size(Pset));
    padj = cell(size(Pset));
    for fam = 1:nFam
        if Ridx(fam)
            [pthr{fam},padj{fam}] = fdrfun(Pset{fam},qval*qfac,cV);
        else
            pthr{fam} = 0;
            padj{fam} = 1;
        end
    end
elseif ismatrix(Pset)
    pthr = zeros(nFam,size(Pset,2));
    padj = zeros(size(Pset));
    for fam = 1:nFam
        if Ridx(fam)
            [pthr(fam),padj(fam,:)] = fdrfun(Pset(fam,:),qval*qfac,cV);
        else
            pthr(fam)   = 0;
            padj(fam,:) = 1;
        end
    end
end

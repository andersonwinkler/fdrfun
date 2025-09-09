function bb2014(varargin)
% Computes the
%
% Note that the corrected and adjusted p-values do **not** depend
% on the supplied q-value, but they do depend on the choice of c(V).
%
% References:
% * Benjamini Y, Bogomolov M. Selective inference on multiple
%   families of hypotheses. Journal of the Royal Statistical
%   Society: Series B (Statistical Methodology).
%   2014 Jan 18;76(1):297–318.
%
% ________________________________
% Anderson M. Winkler
% Research Imaging Center/UTHSCSA
% Dec/2007 (first version)
% Dec/2022 (this version)
% http://brainder.org

narginchk(1,3);
Pset   = varargin{1};
qval   = 0.05;
fdrfun = @bky7;
cV     = 1;
if nargin > 1
    qval   = varargin{2};
end
if nargin > 2
    fdrfun = varargin{2};
end
if nargin == 3
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
        [~,psimes(fam)] = simes(Pset{fam},qval);
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
[~,~,pafos] = fdrfun(psimes,qval,cV); % "p after FDR on Simes"
Ridx = pafos <= qval;
R    = sum(Ridx);

% ========[PART 3: FDR ON EACH SELECTED FAMILY]============================
if iscell(Pset)
    Padj = cell(size(Pset));
    for fam = 1:nFam
        if Ridx(fam)
            [~,~,Padj{fam}] = fdrfun(Pset{fam},qval*R/nFam,cV);
        else
            Padj{fam} = 1;
        end
    end
elseif ismatrix(Pset)
    Padj = zeros(size(Pset));
    for fam = 1:nFam
        if Ridx(fam)
            [~,~,Padj(fam,:)] = fdrfun(Pset(fam,:),qval*R/nFam,cV);
        else
            Padj(fam,:) = 1;
        end
    end
end

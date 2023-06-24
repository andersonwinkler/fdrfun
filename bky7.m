function varargout = bky7(varargin)

% Accept arguments
switch nargin
    case 0
        error('Error: Not enough arguments.');
    case 1
        pval = varargin{1};
        qval = 0.05;
        cV   = 1;
    case 2
        pval = varargin{1};
        qval = varargin{2};
        cV   = 1;
    case 3
        pval = varargin{1};
        qval = varargin{2};
        if varargin{3}
            cV = 1;
        else
            cV = sum(1./(1:numel(pval)));
        end
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

% ========[PART 1: FDR THRESHOLD]========================================

% Sort p-values
[pval,oidx] = sort(pval);

% Number of observations
V = numel(pval);

% Order (indices), in the same size as the pvalues
idx = reshape(1:V,size(pval));

idxthr = false(size(idx));
for v = 1:V
    % Line to be used as cutoff
    thrline = idx(v:end)*qval/(V+1-v*(1-qval))/cV;

    % P-vals that survive the cutoff
    if any(pval(v:end)<=thrline)
        idxthr(v) = true;
    else
        break
    end
end

% Find the largest pval, still under the line
thr = max(pval(idxthr));

% Deal with the case when all the points under the line
% are equal to zero, and other points are above the line
if thr == 0
    v = find(idxthr,1,'last');
    thr = v*qval/(V+1-v*(1-qval))/cV;
end

% Case when it does not cross
if isempty(thr)
    thr = 0;
end

% Returns the result
varargout{1} = thr;

% ========[PART 2: FDR CORRECTED]========================================

if nargout == 2 || nargout == 3
    
    % p-corrected
    pcor = inf(size(pval));
    for v = 1:V
        % For every p-value, this is the minimum q that will eventually
        % satisfy Definition #7 of Benjamini, Krieger and Yekutieli (2006).
        pcor(v) = min(pval(v:end).*(V+1-v).*cV./(idx(v:end)-v*pval(v:end)));
    end

    % Sort back to the original order and output
    [~,oidxR] = sort(oidx);
    varargout{2} = pcor(oidxR);
end

% ========[PART 3: FDR ADJUSTED ]========================================

if nargout == 3

    % The p-adjusted for the current p-value is the cummulative maximum
    % up to it. Note that this is different from equation #3 of
    % Yekutieli & Benjamini (1999) that is used for p-value adjustment of
    % the usual FDR (BH).
    %padj = fliplr(cummin(fliplr(pcor)));
    padj = cummax(pcor);
    varargout{3} = padj(oidxR);
end

% That's it!

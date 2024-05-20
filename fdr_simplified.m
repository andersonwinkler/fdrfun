function [pthr,padj] = fdr_simplified(pval,qval)
[pval,oidx] = sort(pval);
[~,oidxR] = sort(oidx);
V = numel(pval);
idx = reshape(1:V,size(pval));
thrline = idx*qval/V;
pthr = max(pval(pval<=thrline));
if pthr == 0
    pthr = max(thrline(pval<=thrline));
end
if isempty(pthr)
    pthr = 0;
end
pcor = pval.*V./idx;
if isrow(pcor)
    padj = fliplr(cummin(fliplr(pcor)));
else
    padj = flipud(cummin(flipud(pcor)));
end
padj = padj(oidxR);

function [pthr,padj] = myfunc(pval,qval)
[pval,oidx] = sort(pval);
[~,oidxR] = sort(oidx);
V = numel(pval);
idx = reshape(1:V,size(pval));
idxthr = false(size(idx));
for v = 1:V
    thrline = idx(v:end)*qval/(V+1-v*(1-qval));
    if any(pval(v:end)<=thrline)
        idxthr(v) = true;
    else
        break
    end
end
pthr = max(pval(idxthr));
if pthr == 0
    v = find(idxthr,1,'last');
    pthr = v*qval/(V+1-v*(1-qval));
end
if isempty(pthr)
    pthr = 0;
end
pcor = inf(size(pval));
for v = 1:V
    pcor(v) = min(pval(v:end).*(V+1-v)./(idx(v:end)-v*pval(v:end)));
    if pcor(v) >= 1
        pcor(v:end) = 1;
        break
    end
end
padj = cummax(pcor);
padj = padj(oidxR);
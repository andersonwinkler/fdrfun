addpath('~/tools/palm.git')
addpath('../..')
FDRmethod = 'bh1995';

% Choose functions for FDR and for the confidence intervals
switch lower(FDRmethod)
    case {'bh1995','bh'}
        fdrfun = @fdr;
    case {'bky2006','bky'}
        fdrfun = @bky7;
end

% Read the data
T      = palm_miscread('narps-4965_9U7M-hypo1_unthresh.nii',true);
M      = palm_miscread('narps-4965_9U7M-mask.nii',true);
mask   = logical(M.data);
tstats = T.data(mask);
nV     = numel(tstats);

% Indices of positive and negative voxels
idxpos = tstats > 0;
idxneg = ~idxpos;

% The usual pvalues, one-tailed
df           = 107;
pvals        = tcdf(tstats,df,'upper');
T.data(mask) = -log10(pvals);
T.filename   = 'narps-4965_9U7M-hypo1_pvals.nii';
palm_miscwrite(T);
J.pposthr = -log10(.05);
J.pnegthr = 'does not apply';
J.zposthr = -norminv(.05);
J.znegthr = 'does not apply';
writejson(J,'narps-4965_9U7M-hypo1_pvals.json');

% Canonical FDR
[fdrthr,~,fdradj] = fdrfun(pvals);
T.data(mask) = -log10(fdradj);
T.filename   = 'narps-4965_9U7M-hypo1_canonical.nii';
palm_miscwrite(T);
J.pposthr = -log10(fdrthr);
J.pnegthr = 'does not apply';
J.zposthr = -norminv(fdrthr);
J.znegthr = 'does not apply';
writejson(J,'narps-4965_9U7M-hypo1_canonical.json');

% Two-tailed p-values (this is equivalent to what PALM does with -twotail)
pvals2 = 2*tcdf(abs(tstats),df,'upper');
[fdrthr2,~,fdradj2] = fdrfun(pvals2);
T.data(mask)  = -log10(fdradj2);
T.filename    = 'narps-4965_9U7M-hypo1_twotail.nii';
palm_miscwrite(T);
J.pposthr = -log10(fdrthr2);
J.pnegthr = -log10(fdrthr2);
J.zposthr = -norminv(fdrthr2);
J.znegthr = -norminv(fdrthr2);
writejson(J,'narps-4965_9U7M-hypo1_twotail.json');

% Combined pos and neg, with duplicate number of tests (this is what PALM
% does for -corrcon with -fdr, i.e., the cfdrp files)
pvalsc        = [pvals;1-pvals];
[fdrthrc,~,fdradjc] = fdrfun(pvalsc);
T.data(mask)  = -log10(fdradjc(1:nV));
T.filename    = 'narps-4965_9U7M-hypo1_combined_mappos.nii';
palm_miscwrite(T);
T.data(mask)  = -log10(fdradjc(nV+1:end));
T.filename    = 'narps-4965_9U7M-hypo1_combined_mapneg.nii';
palm_miscwrite(T);
J.pposthr = -log10(fdrthrc);
J.pnegthr = -log10(fdrthrc);
J.zposthr = -norminv(fdrthrc);
J.znegthr = -norminv(fdrthrc);
writejson(J,'narps-4965_9U7M-hypo1_combined.json');

% Combined two separate runs of FDR, once all tests, once on all tests negatated
% What Tom had incorrectly had inferred was Chris' suggestion (PALM
% also outputs this even if -corrcon is used, i.e., the fdrp files)
[thr1,~,tmp1] = fdrfun(  pvals);
[thr2,~,tmp2] = fdrfun(1-pvals);
fdradjc2 = [tmp1;tmp2];
T.data(mask)  = -log10(tmp1);
T.filename    = 'narps-4965_9U7M-hypo1_twice_mappos.nii';
palm_miscwrite(T);
T.data(mask)  = -log10(tmp2);
T.filename    = 'narps-4965_9U7M-hypo1_twice_mapneg.nii';
palm_miscwrite(T);
J.pposthr = -log10(thr1);
J.pnegthr = -log10(thr2);
J.zposthr = -norminv(thr1);
J.znegthr = -norminv(thr2);
writejson(J,'narps-4965_9U7M-hypo1_twice.json');

% Combined, but do Sidak first:
psid1 = 1 - (1 - pvals).^2;
psid2 = 1 - (    pvals).^2;
[thr1,~,tmp1] = fdrfun(psid1);
[thr2,~,tmp2] = fdrfun(psid2);
fdradjc3 = [tmp1;tmp2];
T.data(mask)  = -log10(tmp1);
T.filename    = 'narps-4965_9U7M-hypo1_sidaktwice_mappos.nii';
palm_miscwrite(T);
T.data(mask)  = -log10(tmp2);
T.filename    = 'narps-4965_9U7M-hypo1_sidaktwice_mapneg.nii';
palm_miscwrite(T);
J.pposthr = -log10(thr1);
J.pnegthr = -log10(thr2);
J.zposthr = -norminv(thr1);
J.znegthr = -norminv(thr2);
writejson(J,'narps-4965_9U7M-hypo1_sidaktwice.json');

% Only positive or only negative (suggested by Chris)
fdradjs1 = zeros(size(pvals));
[thr1,~,fdradjs1(idxpos)] = fdrfun(  pvals(idxpos));
[thr2,~,fdradjs1(idxneg)] = fdrfun(1-pvals(idxneg));
T.data(mask) = -log10(fdradjs1);
T.filename   = 'narps-4965_9U7M-hypo1_split1tail.nii';
palm_miscwrite(T);
J.pposthr = -log10(thr1);
J.pnegthr = -log10(thr2);
J.zposthr = -norminv(thr1);
J.znegthr = -norminv(thr2);
writejson(J,'narps-4965_9U7M-hypo1_split1tail.json');

% Only positive or only negative (this is Tom's suggestion, but was already
% in Anderson's code as if it had been suggested by Chris)
fdradjs2 = zeros(size(pvals));
[thr1,~,fdradjs2(idxpos)] = fdrfun(pvals2(idxpos));
[thr2,~,fdradjs2(idxneg)] = fdrfun(pvals2(idxneg));
T.data(mask) = -log10(fdradjs2);
T.filename   = 'narps-4965_9U7M-hypo1_split2tail.nii';
palm_miscwrite(T);
J.pposthr = -log10(thr1);
J.pnegthr = -log10(thr2);
J.zposthr = -norminv(thr1);
J.znegthr = -norminv(thr2);
writejson(J,'narps-4965_9U7M-hypo1_split2tail.json');
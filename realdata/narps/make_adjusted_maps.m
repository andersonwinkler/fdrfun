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

% Canonical FDR
[~,~,fdradj] = fdrfun(pvals);
T.data(mask) = -log10(fdradj);
T.filename   = 'narps-4965_9U7M-hypo1_canonical.nii';
palm_miscwrite(T);

% Two-tailed p-values (this is equivalent to what PALM does with -twotail)
pvals2 = 2*tcdf(abs(tstats),df,'upper');
[~,~,fdradj2] = fdrfun(pvals2);
T.data(mask)  = -log10(fdradj2);
T.filename    = 'narps-4965_9U7M-hypo1_twotail.nii';
palm_miscwrite(T);

% Combined pos and neg, with duplicate number of tests (this is what PALM
% does for -corrcon with -fdr, i.e., the cfdrp files)
pvalsc        = [pvals;1-pvals];
[~,~,fdradjc] = fdrfun(pvalsc);
T.data(mask)  = -log10(fdradjc(1:nV));
T.filename    = 'narps-4965_9U7M-hypo1_combined1.nii';
palm_miscwrite(T);
T.data(mask)  = -log10(fdradjc(nV+1:end));
T.filename    = 'narps-4965_9U7M-hypo1_combined2.nii';
palm_miscwrite(T);

% Combined two separate runs of FDR, once all tests, once on all tests negatated
% What Tom had incorrectly had inferred was Chris' suggestion (PALM
% also outputs this even if -corrcon is used, i.e., the fdrp files)
[~,~,tmp1] = fdrfun(  pvals);
[~,~,tmp2] = fdrfun(1-pvals);
fdradjc2 = [tmp1;tmp2];
T.data(mask)  = -log10(tmp1);
T.filename    = 'narps-4965_9U7M-hypo1_twice1.nii';
palm_miscwrite(T);
T.data(mask)  = -log10(tmp2);
T.filename    = 'narps-4965_9U7M-hypo1_twice2.nii';
palm_miscwrite(T);

% Combined, but do Sidak first:
psid1 = 1 - (1 - pvals).^2;
psid2 = 1 - (    pvals).^2;
[~,~,tmp1] = fdrfun(psid1);
[~,~,tmp2] = fdrfun(psid2);
fdradjc3 = [tmp1;tmp2];
T.data(mask)  = -log10(tmp1);
T.filename    = 'narps-4965_9U7M-hypo1_sidaktwice1.nii';
palm_miscwrite(T);
T.data(mask)  = -log10(tmp2);
T.filename    = 'narps-4965_9U7M-hypo1_sidaktwice2.nii';
palm_miscwrite(T);

% Only positive or only negative (suggested by Chris)
fdradjs1 = zeros(size(pvals));
[~,~,fdradjs1(idxpos)] = fdrfun(  pvals(idxpos));
[~,~,fdradjs1(idxneg)] = fdrfun(1-pvals(idxneg));
T.data(mask) = -log10(fdradjs1);
T.filename   = 'narps-4965_9U7M-hypo1_split1tail.nii';
palm_miscwrite(T);

% Only positive or only negative (this is Tom's suggestion, but was already
% in Anderson's code as if it had been suggested by Chris)
fdradjs2 = zeros(size(pvals));
[~,~,fdradjs2(idxpos)] = fdrfun(pvals2(idxpos));
[~,~,fdradjs2(idxneg)] = fdrfun(pvals2(idxneg));
T.data(mask) = -log10(fdradjs2);
T.filename   = 'narps-4965_9U7M-hypo1_split2tail.nii';
palm_miscwrite(T);
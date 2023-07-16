addpath('~/tools/palm.git')
addpath('../..')
FDRmethod = 'bky';

% Choose functions for FDR and for the confidence intervals
switch lower(FDRmethod)
    case {'bh1995','bh'}
        fdrfun = @fdr;
        fdrstr = 'bh';
    case {'bky2006','bky'}
        fdrfun = @bky7;
        fdrstr = 'bky';
end

% Read the data
T      = palm_miscread('stats.FT.nii',true);
M      = palm_miscread('full_mask.FT.nii',true);
mask   = logical(M.data);
tstats = T.data(mask);
nV     = numel(tstats);

% Indices of positive and negative voxels
idxpos = tstats > 0;
idxneg = ~idxpos;

% The usual pvalues, one-tailed
df           = 430;
pvalspos     = tcdf(tstats,df,'upper');
T.data(mask) = -log10(pvalspos);
T.filename   = 'afni_pvals_mappos.nii';
palm_miscwrite(T);
pvalsneg     = tcdf(tstats,df);
T.data(mask) = -log10(pvalsneg);
T.filename   = 'afni_pvals_mapneg.nii';
palm_miscwrite(T);
J.pposthr    = -log10(.05/2);
J.pnegthr    = -log10(.05/2);
J.zposthr    = -tinv(.05/2,df);
J.znegthr    = -tinv(.05/2,df);
writejson(J,'afni_pvals.json');

% Canonical FDR
[fdrthr,~,fdradj] = fdrfun(pvalspos);
T.data(mask) = -log10(fdradj);
T.filename   = sprintf('afni_%s_canonical.nii',fdrstr);
palm_miscwrite(T);
J.pposthr    = -log10(fdrthr);
J.pnegthr    = 'does not apply';
J.zposthr    = -tinv(fdrthr,df);
J.znegthr    = 'does not apply';
writejson(J,sprintf('afni_%s_canonical.json',fdrstr));

% Two-tailed p-values (this is equivalent to what PALM does with -twotail)
pvals2 = 2*tcdf(abs(tstats),df,'upper');
[fdrthr2,~,fdradj2] = fdrfun(pvals2);
T.data(mask)  = -log10(fdradj2);
T.filename    = sprintf('afni_%s_twotail.nii',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(fdrthr2);
J.pnegthr     = -log10(fdrthr2);
J.zposthr     = -tinv(fdrthr2,df);
J.znegthr     = -tinv(fdrthr2,df);
writejson(J,sprintf('afni_%s_twotail.json',fdrstr));

% Combined pos and neg, with duplicate number of tests (this is what PALM
% does for -corrcon with -fdr, i.e., the cfdrp files)
pvalsc        = [pvalspos;pvalsneg];
[fdrthrc,~,fdradjc] = fdrfun(pvalsc);
T.data(mask)  = -log10(fdradjc(1:nV));
T.filename    = sprintf('afni_%s_combined_mappos.nii',fdrstr);
palm_miscwrite(T);
T.data(mask)  = -log10(fdradjc(nV+1:end));
T.filename    = sprintf('afni_%s_combined_mapneg.nii',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(fdrthrc);
J.pnegthr     = -log10(fdrthrc);
J.zposthr     = -tinv(fdrthrc,df);
J.znegthr     = -tinv(fdrthrc,df);
writejson(J,sprintf('afni_%s_combined.json',fdrstr));

% Combined two separate runs of FDR, once all tests, once on all tests negatated
% What Tom had incorrectly had inferred was Chris' suggestion (PALM
% also outputs this even if -corrcon is used, i.e., the fdrp files)
[thr1,~,tmp1] = fdrfun(pvalspos);
[thr2,~,tmp2] = fdrfun(pvalsneg);
fdradjc2 = [tmp1;tmp2];
T.data(mask)  = -log10(tmp1);
T.filename    = sprintf('afni_%s_twice_mappos.nii',fdrstr);
palm_miscwrite(T);
T.data(mask)  = -log10(tmp2);
T.filename    = sprintf('afni_%s_twice_mapneg.nii',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thr1);
J.pnegthr     = -log10(thr2);
J.zposthr     = -tinv(thr1,df);
J.znegthr     = -tinv(thr2,df);
writejson(J,sprintf('afni_%s_twice.json',fdrstr));

% Combined, but do Sidak first:
psid1 = 1 - (1 - pvalspos).^2;
psid2 = 1 - (1 - pvalsneg).^2;
[thr1,~,tmp1] = fdrfun(psid1);
[thr2,~,tmp2] = fdrfun(psid2);
fdradjc3 = [tmp1;tmp2];
T.data(mask)  = -log10(tmp1);
T.filename    = sprintf('afni_%s_sidaktwice_mappos.nii',fdrstr);
palm_miscwrite(T);
T.data(mask)  = -log10(tmp2);
T.filename    = sprintf('afni_%s_sidaktwice_mapneg.nii',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thr1);
J.pnegthr     = -log10(thr2);
J.zposthr     = -tinv(thr1,df);
J.znegthr     = -tinv(thr2,df);
writejson(J,sprintf('afni_%s_sidaktwice.json',fdrstr));

% Only positive or only negative (suggested by Chris)
fdradjs1 = zeros(size(pvalspos));
[thr1,~,fdradjs1(idxpos)] = fdrfun(pvalspos(idxpos));
[thr2,~,fdradjs1(idxneg)] = fdrfun(pvalsneg(idxneg));
T.data(mask) = -log10(fdradjs1);
T.filename   = sprintf('afni_%s_split1tail.nii',fdrstr);
palm_miscwrite(T);
J.pposthr    = -log10(thr1);
J.pnegthr    = -log10(thr2);
J.zposthr    = -tinv(thr1,df);
J.znegthr    = -tinv(thr2,df);
writejson(J,sprintf('afni_%s_split1tail.json',fdrstr));

% Only positive or only negative (this is Tom's suggestion, but was already
% in Anderson's code as if it had been suggested by Chris)
fdradjs2 = zeros(size(pvalspos));
[thr1,~,fdradjs2(idxpos)] = fdrfun(pvals2(idxpos));
[thr2,~,fdradjs2(idxneg)] = fdrfun(pvals2(idxneg));
T.data(mask) = -log10(fdradjs2);
T.filename   = sprintf('afni_%s_split2tail.nii',fdrstr);
palm_miscwrite(T);
J.pposthr    = -log10(thr1);
J.pnegthr    = -log10(thr2);
J.zposthr    = -tinv(thr1,df);
J.znegthr    = -tinv(thr2,df);
writejson(J,sprintf('afni_%s_split2tail.json',fdrstr));
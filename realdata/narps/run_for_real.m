addpath('~/tools/palm.git')
addpath('../..')
FDRmethod = 'bh';
alpha = 0.05; % our q-threshold

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
T      = palm_miscread('narps-4965_9U7M-hypo1_unthresh.nii',false);
M      = palm_miscread('narps-4965_9U7M-mask.nii',false);
mask   = logical(M.data);
tstats = T.data(mask);
nV     = numel(tstats);

% Indices of positive and negative voxels
idxpos = tstats > 0;
idxneg = ~idxpos;

% UNCORRECTED
df           = 107;
pvals0       = tcdf(tstats,df,'upper');
T.data(mask) = -log10(pvals0);
T.filename   = 'narps-4965_9U7M-hypo1_unc_pvals_mappos';
palm_miscwrite(T);
pvals1       = tcdf(tstats,df); % 1-p, to avoid losing precision
T.data(mask) = -log10(pvals1);
T.filename   = 'narps-4965_9U7M-hypo1_unc_pvals_mapneg';
palm_miscwrite(T);
J.pposthr    = -log10(alpha/2);
J.pnegthr    = -log10(alpha/2);
J.zposthr    = -tinv (alpha/2,df);
J.znegthr    = -tinv (alpha/2,df);
writejson(J,'narps-4965_9U7M-hypo1_unc_pvals.json');

% CANONICAL
[thrcan0,adjcan0] = fdrfun(pvals0);
[thrcan1,adjcan1] = fdrfun(pvals1);
T.data(mask)  = -log10(adjcan0);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_canonical_mappos',fdrstr);
palm_miscwrite(T);
T.data(mask)  = -log10(adjcan1);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_canonical_mapneg',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thrcan0);
J.pnegthr     = -log10(thrcan1);
J.zposthr     = -tinv( thrcan0,df);
J.znegthr     = -tinv( thrcan1,df);
writejson(J,sprintf('narps-4965_9U7M-hypo1_%s_canonical.json',fdrstr));

% COMBINED
[thrcom,adjcom] = fdrfun([pvals0;pvals1]);
T.data(mask)  = -log10(adjcom(1:nV));
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_combined_mappos',fdrstr);
palm_miscwrite(T);
T.data(mask)  = -log10(adjcom(nV+1:end));
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_combined_mapneg',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thrcom);
J.pnegthr     = -log10(thrcom);
J.zposthr     = -tinv( thrcom,df);
J.znegthr     = -tinv( thrcom,df);
writejson(J,sprintf('narps-4965_9U7M-hypo1_%s_combined.json',fdrstr));

% TWO-TAILED
pvals2 = 2*tcdf(abs(tstats),df,'upper');
[thrtwo,adjtwo] = fdrfun(pvals2);
T.data(mask)  = -log10(adjtwo);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_twotail',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thrtwo);
J.pnegthr     = -log10(thrtwo);
J.zposthr     = -tinv( thrtwo,df);
J.znegthr     = -tinv( thrtwo,df);
writejson(J,sprintf('narps-4965_9U7M-hypo1_%s_twotail.json',fdrstr));

% SPLIT + TWO-TAILED
adjspl = ones(size(pvals0));
[thrspl0,adjspl(idxpos)] = fdrfun(pvals2(idxpos));
[thrspl1,adjspl(idxneg)] = fdrfun(pvals2(idxneg));
T.data(mask)  = -log10(adjspl);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_split2tail',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thrspl0);
J.pnegthr     = -log10(thrspl1);
J.zposthr     = -tinv( thrspl0,df);
J.znegthr     = -tinv( thrspl1,df);
writejson(J,sprintf('narps-4965_9U7M-hypo1_%s_split2tail.json',fdrstr));

% CANONICAL + BB2014
Pset = [pvals0';pvals1'];
[thrcanbb,adjcanbb] = bb2014(Pset,[],fdrfun);
thrcanbb0 = thrcanbb(1);
thrcanbb1 = thrcanbb(2);
adjcanbb0 = adjcanbb(1,:)';
adjcanbb1 = adjcanbb(2,:)';
T.data(mask)  = -log10(adjcanbb0);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_canonicalbb_mappos',fdrstr);
palm_miscwrite(T);
T.data(mask)  = -log10(adjcanbb1);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_canonicalbb_mapneg',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thrcanbb0);
J.pnegthr     = -log10(thrcanbb1);
J.zposthr     = -tinv( thrcanbb0,df);
J.znegthr     = -tinv( thrcanbb1,df);
writejson(J,sprintf('narps-4965_9U7M-hypo1_%s_canonicalbb.json',fdrstr));

% SPLIT + TWO-TAILED + BB2014
adjsplbb = ones(size(pvals0));
Pset = {pvals2(idxpos), pvals2(idxneg)};
[tmpthr,tmpadj] = bb2014(Pset,[],fdrfun);
thrsplbb0 = tmpthr{1};
thrsplbb1 = tmpthr{2};
adjsplbb(idxpos) = tmpadj{1};
adjsplbb(idxneg) = tmpadj{2};
T.data(mask)  = -log10(adjsplbb);
T.filename    = sprintf('narps-4965_9U7M-hypo1_%s_split2tailbb',fdrstr);
palm_miscwrite(T);
J.pposthr     = -log10(thrsplbb0);
J.pnegthr     = -log10(thrsplbb1);
J.zposthr     = -tinv( thrsplbb0,df);
J.znegthr     = -tinv( thrsplbb1,df);
writejson(J,sprintf('narps-4965_9U7M-hypo1_%s_split2tailbb.json',fdrstr));
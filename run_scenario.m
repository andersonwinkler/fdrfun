function run_scenario(commonfile,scenariofile,outputfile)

% Read general configuration (common to all scenarios)
J = readjson(commonfile);

% Read json with scenario configuration
Js = readjson(scenariofile);

% Merge them, for the outputs later
F = fieldnames(Js);
for f = 1:numel(F)
    J.(F{f}) = Js.(F{f});
end

% Simulation parameters
if isoctave()
    rand('seed',J.seed) %#ok<RAND>
else
    rng(J.seed);
end
numTests        = J.numTests; % number of "voxels"
effectSize      = J.effectSize; % z-stat that determines the synthetic effect size
fracPos         = J.fracPos; % fraction of tests with true positive effect
fracNeg         = J.fracNeg; % fraction of tests with true negative effect
rho             = J.rho;  % compound symmetric correlation among the numTests (rho>=0)
q               = J.q; % test level, E(FDR) to be controlled
numRealizations = J.numRealizations; % number of times we repeat the simulation
FDRmethod       = J.FDRmethod; % use 'bh1995' or 'bky2006'
CImethod        = J.CImethod; % use 'Wald' or 'Wilson'
alpha           = J.alpha; % for the confidence interval

% Vars for later
fdp_fdr       = zeros(numRealizations,1);
fdp_fdr2      = zeros(numRealizations,1);
fdp_fdrc      = zeros(numRealizations,1);
fdp_fdrc2     = zeros(numRealizations,1);
fdp_fdrc3     = zeros(numRealizations,1);
fdp_fdrs1     = zeros(numRealizations,1);
fdp_fdrs2     = zeros(numRealizations,1);
fdp_fdr_pos   = zeros(numRealizations,1);
fdp_fdr2_pos  = zeros(numRealizations,1);
fdp_fdrc_pos  = zeros(numRealizations,1);
fdp_fdrc2_pos = zeros(numRealizations,1);
fdp_fdrc3_pos = zeros(numRealizations,1);
fdp_fdrs1_pos = zeros(numRealizations,1);
fdp_fdrs2_pos = zeros(numRealizations,1);
fdp_fdr_neg   = zeros(numRealizations,1);
fdp_fdr2_neg  = zeros(numRealizations,1);
fdp_fdrc_neg  = zeros(numRealizations,1);
fdp_fdrc2_neg = zeros(numRealizations,1);
fdp_fdrc3_neg = zeros(numRealizations,1);
fdp_fdrs1_neg = zeros(numRealizations,1);
fdp_fdrs2_neg = zeros(numRealizations,1);

% Choose functions for FDR and for the confidence intervals
switch lower(FDRmethod)
    case {'bh1995','bh'}
        fdrfun = @fdr;
    case {'bky2006','bky'}
        fdrfun = @bky7;
end
switch lower(CImethod)
    case 'wald'
        confint = @(x) mean(x,1) + [-1,1]*norminv(1-alpha/2)*std(x,[],1)/sqrt(size(x,1));
    case 'wilson'
        confint = @(x) wilson(size(x,1),sum(x,1),alpha);
end

% Signal to be added and where
numPos  = round(numTests * fracPos);
numNeg  = round(numTests * fracNeg);
maskPos = (1:numTests)'    <= numPos;
maskNeg = (numTests:-1:1)' <= numNeg;
signal  = zeros(numTests,1);
if numPos >= 1
    signal(1:numPos) = effectSize;
end
if numNeg >= 1
    signal(end-numNeg:end) = -effectSize;
end

% For each realization
for rlz = 1:numRealizations

    % Create random data, add signal
    zstats = signal + sqrt(1-rho)*randn(numTests,1) + sqrt(rho)*randn(1,1);
    idxpos = zstats > 0;
    idxneg = ~idxpos;

    % The usual pvalues, one-tailed
    pvals = normcdf(zstats,'upper');
    [~,~,fdradj] = fdrfun(pvals);

    % Two-tailed p-values (this is equivalent to what PALM does with -twotail)
    pvals2 = 2*normcdf(abs(zstats),'upper');
    [~,~,fdradj2] = fdrfun(pvals2);

    % Combined pos and neg, with duplicate number of tests (this is what PALM
    % does for -corrcon with -fdr, i.e., the cfdrp files)
    pvalsc = [pvals;1-pvals];
    [~,~,fdradjc] = fdrfun(pvalsc);

    % Combined two separate runs of FDR, once all tests, once on all tests negatated
    % What Tom had incorrectly had inferred was Chris' suggestion (PALM
    % also outputs this even if -corrcon is used, i.e., the fdrp files)
    [~,~,tmp1] = fdrfun(  pvals);
    [~,~,tmp2] = fdrfun(1-pvals);
    fdradjc2 = [tmp1;tmp2];

    % Combined, but do Sidak first:
    psid1 = 1 - (1 - pvals).^2;
    psid2 = 1 - (    pvals).^2;
    [~,~,tmp1] = fdrfun(psid1);
    [~,~,tmp2] = fdrfun(psid2);
    fdradjc3 = [tmp1;tmp2];

    % Only positive or only negative (suggested by Chris)
    fdradjs1 = zeros(size(pvals));
    [~,~,fdradjs1(idxpos)] = fdrfun(  pvals(idxpos));
    [~,~,fdradjs1(idxneg)] = fdrfun(1-pvals(idxneg));

    % Only positive or only negative (this is Tom's suggestion, but was already
    % in Anderson's code as if it had been suggested by Chris)
    fdradjs2 = zeros(size(pvals));
    [~,~,fdradjs2(idxpos)] = fdrfun(pvals2(idxpos));
    [~,~,fdradjs2(idxneg)] = fdrfun(pvals2(idxneg));

    % Empirical FDRs (global, i.e., the user looks into both sides of the map)
    fdp_fdr   (rlz) = sum((fdradj   <= q) & ~  maskPos)          / sum(fdradj   <= q);
    fdp_fdr2  (rlz) = sum((fdradj2  <= q) & ~ (maskPos|maskNeg)) / sum(fdradj2  <= q);
    fdp_fdrc  (rlz) = sum((fdradjc  <= q) & ~ [maskPos;maskNeg]) / sum(fdradjc  <= q);
    fdp_fdrc2 (rlz) = sum((fdradjc2 <= q) & ~ [maskPos;maskNeg]) / sum(fdradjc2 <= q);
    fdp_fdrc3 (rlz) = sum((fdradjc3 <= q) & ~ [maskPos;maskNeg]) / sum(fdradjc3 <= q);
    fdp_fdrs1 (rlz) = sum((fdradjs1 <= q) & ~ (maskPos|maskNeg)) / sum(fdradjs1 <= q);
    fdp_fdrs2 (rlz) = sum((fdradjs2 <= q) & ~ (maskPos|maskNeg)) / sum(fdradjs2 <= q);

    % Empirical FDRs (positive side of the map, i.e., user interested in positive results only)
    idxpoc = [idxpos;false(size(idxpos))];
    fdp_fdr_pos   (rlz) = sum((fdradj   (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradj   (idxpos) <= q);
    fdp_fdr2_pos  (rlz) = sum((fdradj2  (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradj2  (idxpos) <= q);
    fdp_fdrc_pos  (rlz) = sum((fdradjc  (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(fdradjc  (idxpoc) <= q);
    fdp_fdrc2_pos (rlz) = sum((fdradjc2 (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(fdradjc2 (idxpoc) <= q);
    fdp_fdrc3_pos (rlz) = sum((fdradjc3 (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(fdradjc3 (idxpoc) <= q);
    fdp_fdrs1_pos (rlz) = sum((fdradjs1 (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradjs1 (idxpos) <= q);
    fdp_fdrs2_pos (rlz) = sum((fdradjs2 (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradjs2 (idxpos) <= q);

    % Empirical FDRs (negative side of the map, i.e., user interested in negative results only)
    idxnec = [false(size(idxneg));idxneg];
    fdp_fdr_neg   (rlz) = sum((fdradj   (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradj   (idxneg) <= q);
    fdp_fdr2_neg  (rlz) = sum((fdradj2  (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradj2  (idxneg) <= q);
    fdp_fdrc_neg  (rlz) = sum((fdradjc  (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(fdradjc  (idxnec) <= q);
    fdp_fdrc2_neg (rlz) = sum((fdradjc2 (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(fdradjc2 (idxnec) <= q);
    fdp_fdrc3_neg (rlz) = sum((fdradjc3 (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(fdradjc3 (idxnec) <= q);
    fdp_fdrs1_neg (rlz) = sum((fdradjs1 (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradjs1 (idxneg) <= q);
    fdp_fdrs2_neg (rlz) = sum((fdradjs2 (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradjs2 (idxneg) <= q);
end

% If the FDP is NaN, then there were no positives, which we'd interpret as
% zero false positives (as opposed to 0/0 false positives)
fdp_fdr       (isnan(fdp_fdr  )) = 0;
fdp_fdr2      (isnan(fdp_fdr2 )) = 0;
fdp_fdrc      (isnan(fdp_fdrc )) = 0;
fdp_fdrc2     (isnan(fdp_fdrc2)) = 0;
fdp_fdrc3     (isnan(fdp_fdrc3)) = 0;
fdp_fdrs1     (isnan(fdp_fdrs1)) = 0;
fdp_fdrs2     (isnan(fdp_fdrs2)) = 0;
fdp_fdr_pos   (isnan(fdp_fdr_pos  )) = 0;
fdp_fdr2_pos  (isnan(fdp_fdr2_pos )) = 0;
fdp_fdrc_pos  (isnan(fdp_fdrc_pos )) = 0;
fdp_fdrc2_pos (isnan(fdp_fdrc2_pos)) = 0;
fdp_fdrc3_pos (isnan(fdp_fdrc3_pos)) = 0;
fdp_fdrs1_pos (isnan(fdp_fdrs1_pos)) = 0;
fdp_fdrs2_pos (isnan(fdp_fdrs2_pos)) = 0;
fdp_fdr_neg   (isnan(fdp_fdr_neg  )) = 0;
fdp_fdr2_neg  (isnan(fdp_fdr2_neg )) = 0;
fdp_fdrc_neg  (isnan(fdp_fdrc_neg )) = 0;
fdp_fdrc2_neg (isnan(fdp_fdrc2_neg)) = 0;
fdp_fdrc3_neg (isnan(fdp_fdrc3_neg)) = 0;
fdp_fdrs1_neg (isnan(fdp_fdrs1_neg)) = 0;
fdp_fdrs2_neg (isnan(fdp_fdrs2_neg)) = 0;

% Compute results
J.BothSides.canonical  = [mean(fdp_fdr),       confint(fdp_fdr)];
J.BothSides.twotailed  = [mean(fdp_fdr2),      confint(fdp_fdr2)];
J.BothSides.combined   = [mean(fdp_fdrc),      confint(fdp_fdrc)];
J.BothSides.twice      = [mean(fdp_fdrc2),     confint(fdp_fdrc2)];
J.BothSides.sidaktwice = [mean(fdp_fdrc3),     confint(fdp_fdrc3)];
J.BothSides.split1tail = [mean(fdp_fdrs1),     confint(fdp_fdrs1)];
J.BothSides.split2tail = [mean(fdp_fdrs2),     confint(fdp_fdrs2)];
J.PosSide.canonical    = [mean(fdp_fdr_pos),   confint(fdp_fdr_pos)];
J.PosSide.twotailed    = [mean(fdp_fdr2_pos),  confint(fdp_fdr2_pos)];
J.PosSide.combined     = [mean(fdp_fdrc_pos),  confint(fdp_fdrc_pos)];
J.PosSide.twice        = [mean(fdp_fdrc2_pos), confint(fdp_fdrc2_pos)];
J.PosSide.sidaktwice   = [mean(fdp_fdrc3_pos), confint(fdp_fdrc3_pos)];
J.PosSide.split1tail   = [mean(fdp_fdrs1_pos), confint(fdp_fdrs1_pos)];
J.PosSide.split2tail   = [mean(fdp_fdrs2_pos), confint(fdp_fdrs2_pos)];
J.NegSide.canonical    = [mean(fdp_fdr_neg),   confint(fdp_fdr_neg)];
J.NegSide.twotailed    = [mean(fdp_fdr2_neg),  confint(fdp_fdr2_neg)];
J.NegSide.combined     = [mean(fdp_fdrc_neg),  confint(fdp_fdrc_neg)];
J.NegSide.twice        = [mean(fdp_fdrc2_neg), confint(fdp_fdrc2_neg)];
J.NegSide.sidaktwice   = [mean(fdp_fdrc3_neg), confint(fdp_fdrc3_neg)];
J.NegSide.split1tail   = [mean(fdp_fdrs1_neg), confint(fdp_fdrs1_neg)];
J.NegSide.split2tail   = [mean(fdp_fdrs2_neg), confint(fdp_fdrs2_neg)];

% Write results
writejson(J,outputfile)

% Simulation parameters
rng('shuffle');
numTests   = 1000;
effectSize = 3;
fracPos    = 0.1;
fracNeg    = 0.4;
rho        = 0.8;  % Compound symmetric correlation among the numTests (rho>=0)
q          = 0.05; % test level, E(FDR) to be controlled
nRlz       = 1000;
fdr_method = 'bky2006'; % use 'bh1995' or 'bky2006'
ci_method  = 'Wald'; % use 'Wald' or 'Wilson'
alpha      = 0.05; % for the confidence interval

% Vars for later
fdp_fdr   = zeros(nRlz,1);
fdp_fdr2  = zeros(nRlz,1);
fdp_fdrc  = zeros(nRlz,1);
fdp_fdrc2 = zeros(nRlz,1);
fdp_fdrs  = zeros(nRlz,1);
fdp_fdrs2 = zeros(nRlz,1);
fdp_fdr_pos   = zeros(nRlz,1);
fdp_fdr2_pos  = zeros(nRlz,1);
fdp_fdrc_pos  = zeros(nRlz,1);
fdp_fdrc2_pos = zeros(nRlz,1);
fdp_fdrs_pos  = zeros(nRlz,1);
fdp_fdrs2_pos = zeros(nRlz,1);
fdp_fdr_neg   = zeros(nRlz,1);
fdp_fdr2_neg  = zeros(nRlz,1);
fdp_fdrc_neg  = zeros(nRlz,1);
fdp_fdrc2_neg = zeros(nRlz,1);
fdp_fdrs_neg  = zeros(nRlz,1);
fdp_fdrs2_neg = zeros(nRlz,1);

% Choose functions for FDR and for the confidence intervals
switch lower(fdr_method)
    case {'bh1995','bh'}
        fdrfun = @fdr;
    case {'bky2006','bky'}
        fdrfun = @bky7;
end
switch lower(ci_method)
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
for rlz = 1:nRlz

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

    % Only positive or only negative (suggested by Chris)
    fdradjs = zeros(size(pvals));
    [~,~,fdradjs(idxpos)] = fdrfun(  pvals(idxpos));
    [~,~,fdradjs(idxneg)] = fdrfun(1-pvals(idxneg));

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
    fdp_fdrs  (rlz) = sum((fdradjs  <= q) & ~ (maskPos|maskNeg)) / sum(fdradjs  <= q);
    fdp_fdrs2 (rlz) = sum((fdradjs2 <= q) & ~ (maskPos|maskNeg)) / sum(fdradjs2 <= q);

    % Empirical FDRs (positive side of the map, i.e., user interested in positive results only)
    idxpoc = [idxpos;false(size(idxpos))];
    fdp_fdr_pos   (rlz) = sum((fdradj   (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradj   (idxpos) <= q);
    fdp_fdr2_pos  (rlz) = sum((fdradj2  (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradj2  (idxpos) <= q);
    fdp_fdrc_pos  (rlz) = sum((fdradjc  (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(fdradjc  (idxpoc) <= q);
    fdp_fdrc2_pos (rlz) = sum((fdradjc2 (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(fdradjc2 (idxpoc) <= q);
    fdp_fdrs_pos  (rlz) = sum((fdradjs  (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradjs  (idxpos) <= q);
    fdp_fdrs2_pos (rlz) = sum((fdradjs2 (idxpos) <= q) & ~ maskPos(idxpos)) / sum(fdradjs2 (idxpos) <= q);

    % Empirical FDRs (negative side of the map, i.e., user interested in negative results only)
    idxnec = [false(size(idxneg));idxneg];
    fdp_fdr_neg   (rlz) = sum((fdradj   (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradj   (idxneg) <= q);
    fdp_fdr2_neg  (rlz) = sum((fdradj2  (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradj2  (idxneg) <= q);
    fdp_fdrc_neg  (rlz) = sum((fdradjc  (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(fdradjc  (idxnec) <= q);
    fdp_fdrc2_neg (rlz) = sum((fdradjc2 (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(fdradjc2 (idxnec) <= q);
    fdp_fdrs_neg  (rlz) = sum((fdradjs  (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradjs  (idxneg) <= q);
    fdp_fdrs2_neg (rlz) = sum((fdradjs2 (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(fdradjs2 (idxneg) <= q);
end

% If the FDP is NaN, then there were no positives, which we'd interpret as
% zero false positives (as opposed to 0/0 false positives)
fdp_fdr  (isnan(fdp_fdr  )) = 0;
fdp_fdr2 (isnan(fdp_fdr2 )) = 0;
fdp_fdrc (isnan(fdp_fdrc )) = 0;
fdp_fdrc2(isnan(fdp_fdrc2)) = 0;
fdp_fdrs (isnan(fdp_fdrs )) = 0;
fdp_fdrs2(isnan(fdp_fdrs2)) = 0;
fdp_fdr_pos  (isnan(fdp_fdr_pos  )) = 0;
fdp_fdr2_pos (isnan(fdp_fdr2_pos )) = 0;
fdp_fdrc_pos (isnan(fdp_fdrc_pos )) = 0;
fdp_fdrc2_pos(isnan(fdp_fdrc2_pos)) = 0;
fdp_fdrs_pos (isnan(fdp_fdrs_pos )) = 0;
fdp_fdrs2_pos(isnan(fdp_fdrs2_pos)) = 0;
fdp_fdr_neg  (isnan(fdp_fdr_neg  )) = 0;
fdp_fdr2_neg (isnan(fdp_fdr2_neg )) = 0;
fdp_fdrc_neg (isnan(fdp_fdrc_neg )) = 0;
fdp_fdrc2_neg(isnan(fdp_fdrc2_neg)) = 0;
fdp_fdrs_neg (isnan(fdp_fdrs_neg )) = 0;
fdp_fdrs2_neg(isnan(fdp_fdrs2_neg)) = 0;

% Print a report
fprintf('Simulation parameters:\n');
fprintf('- Number of tests: %d\n',numTests);
fprintf('- Effect size (z): %g\n',effectSize);
fprintf('- Fraction of positive tests: %g\n',fracPos);
fprintf('- Fraction of negative tests: %g\n',fracNeg);
fprintf('- Compound symmetric correlation: %g\n',rho);
fprintf('- FDR method: %s\n',upper(fdr_method));
fprintf('- Number of realizations: %d\n',nRlz);
fprintf('- Confidence interval: %g%% / Method: %s\n',100*(1-alpha),ci_method);
fprintf('\n');
fprintf('Side of the map: BOTH\n');
fprintf('- E(FDP) canonical FDR:      %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdr),  confint(fdp_fdr));
fprintf('- E(FDP) two-tailed FDR:     %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdr2), confint(fdp_fdr2));
fprintf('- E(FDP) combined FDR:       %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrc), confint(fdp_fdrc));
fprintf('- E(FDP) FDR twice:          %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrc2),confint(fdp_fdrc2));
fprintf('- E(FDP) split + 1tail FDR:  %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrs), confint(fdp_fdrs));
fprintf('- E(FDP) split + 2tail FDR:  %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrs2),confint(fdp_fdrs2));
fprintf('\n');
fprintf('Side of the map: POSITIVE\n');
fprintf('- E(FDP) canonical FDR:      %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdr_pos),  confint(fdp_fdr_pos));
fprintf('- E(FDP) two-tailed FDR:     %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdr2_pos), confint(fdp_fdr2_pos));
fprintf('- E(FDP) combined FDR:       %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrc_pos), confint(fdp_fdrc_pos));
fprintf('- E(FDP) FDR twice:          %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrc2_pos),confint(fdp_fdrc2_pos));
fprintf('- E(FDP) split + 1tail FDR:  %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrs_pos), confint(fdp_fdrs_pos));
fprintf('- E(FDP) split + 2tail FDR:  %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrs2_pos),confint(fdp_fdrs2_pos));
fprintf('\n');
fprintf('Side of the map: NEGATIVE\n');
fprintf('- E(FDP) canonical FDR:      %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdr_neg),  confint(fdp_fdr_neg));
fprintf('- E(FDP) two-tailed FDR:     %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdr2_neg), confint(fdp_fdr2_neg));
fprintf('- E(FDP) combined FDR:       %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrc_neg), confint(fdp_fdrc_neg));
fprintf('- E(FDP) FDR twice:          %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrc2_neg),confint(fdp_fdrc2_neg));
fprintf('- E(FDP) split + 1tail FDR:  %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrs_neg), confint(fdp_fdrs_neg));
fprintf('- E(FDP) split + 2tail FDR:  %0.4f (%0.4f-%0.4f)\n',mean(fdp_fdrs2_neg),confint(fdp_fdrs2_neg));
fprintf('\n');

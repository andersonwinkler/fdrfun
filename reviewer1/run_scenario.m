function run_scenario(commonfile,scenariofile,outputfile)

% Keep track of time
st = tic;

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
numTests        = J.numTests;        % number of "voxels"
effectSize      = J.effectSize;      % z-stat that determines the synthetic effect size
fracPos         = J.fracPos;         % fraction of tests with true positive effect
fracNeg         = J.fracNeg;         % fraction of tests with true negative effect
rho_sign        = J.rho_sign;        % sign of the compound symmetric correlation among the numTests
q               = J.q;               % test level, E(FDR) to be controlled
numRealizations = J.numRealizations; % number of times we repeat the simulation
FDRmethod       = J.FDRmethod;       % use 'bh1995' or 'bky2006'
CImethod        = J.CImethod;        % use 'Wald' or 'Wilson'
alpha           = J.alpha;           % for the confidence interval

% Vars for later
fdp_can     = zeros(numRealizations,1);
fdp_com     = zeros(numRealizations,1);
fdp_two     = zeros(numRealizations,1);
fdp_spl     = zeros(numRealizations,1);
fdp_Old     = zeros(numRealizations,1);
fdp_Q11     = zeros(numRealizations,1);
fdp_Q12     = zeros(numRealizations,1);
fdp_can_pos = zeros(numRealizations,1);
fdp_com_pos = zeros(numRealizations,1);
fdp_two_pos = zeros(numRealizations,1);
fdp_spl_pos = zeros(numRealizations,1);
fdp_Old_pos = zeros(numRealizations,1);
fdp_Q11_pos = zeros(numRealizations,1);
fdp_Q12_pos = zeros(numRealizations,1);
fdp_can_neg = zeros(numRealizations,1);
fdp_com_neg = zeros(numRealizations,1);
fdp_two_neg = zeros(numRealizations,1);
fdp_spl_neg = zeros(numRealizations,1);
fdp_Old_neg = zeros(numRealizations,1);
fdp_Q11_neg = zeros(numRealizations,1);
fdp_Q12_neg = zeros(numRealizations,1);

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

% Compound symmetric correlation among the numTests
rho = rho_sign/numTests;
C   = ones(numTests)*rho + diag(ones(numTests,1)*(1-rho));
R   = chol(C);

% For each realization
for rlz = 1:numRealizations

    % Create random data, add signal
    zstats = signal + R'*randn(numTests,1);
    idxpos = zstats > 0;
    idxneg = ~idxpos;

    % The usual pvalues, one-tailed
    pvals0 = normcdf(zstats,'upper');
    pvals1 = normcdf(zstats); % This is to avoid computing 1-p later and thus retain precision

    % Two-tailed p-values
    pvals2 = 2*normcdf(abs(zstats),'upper');

    % CANONICAL
    [~,~,adj1] = fdrfun(pvals0);
    [~,~,adj2] = fdrfun(pvals1);
    adjcan = [adj1;adj2];

    % COMBINED
    [~,~,adjcom] = fdrfun([pvals0;pvals1]);

    % TWO-TAILED
    [~,~,adjtwo] = fdrfun(pvals2);

    % SPLIT + TWO-TAILED
    adjspl = zeros(size(pvals0));
    [~,~,adjspl(idxpos)] = fdrfun(pvals2(idxpos));
    [~,~,adjspl(idxneg)] = fdrfun(pvals2(idxneg));

    % FOR REVIEWER 1 - OLDER RESULT NOT SHOWN IN THE MANUSCRIPT
    % Similar to combined, but do Sidak first (could have been Bonferroni)
    psid1 = 1 - (pvals1).^2;
    psid0 = 1 - (pvals0).^2;
    [~,~,tmp1] = fdrfun(psid1);
    [~,~,tmp2] = fdrfun(psid0);
    adjOld = [tmp1;tmp2];

    % REVIEWER 1 - METHOD A - Question 11
    adjQ11 = adjspl*2;

    % REVIEWER 1 - METHOD B - Question 12
    adjQ12 = adjcan*2;

    % Empirical FDRs (global, i.e., the user looks into both sides of the map)
    fdp_can(rlz) = sum((adjcan <= q) & ~ [maskPos;maskNeg]) / sum(adjcan <= q);
    fdp_com(rlz) = sum((adjcom <= q) & ~ [maskPos;maskNeg]) / sum(adjcom <= q);
    fdp_two(rlz) = sum((adjtwo <= q) & ~ (maskPos|maskNeg)) / sum(adjtwo <= q);
    fdp_spl(rlz) = sum((adjspl <= q) & ~ (maskPos|maskNeg)) / sum(adjspl <= q);
    fdp_Old(rlz) = sum((adjOld <= q) & ~ [maskPos;maskNeg]) / sum(adjOld <= q);
    fdp_Q11(rlz) = sum((adjQ11 <= q) & ~ (maskPos|maskNeg)) / sum(adjQ11 <= q);
    fdp_Q12(rlz) = sum((adjQ12 <= q) & ~ [maskPos;maskNeg]) / sum(adjQ12 <= q);

    % Empirical FDRs (positive side of the map, i.e., user interested in positive results only, i.e., results that match the direction of the contrast)
    idxpoc = [idxpos;false(size(idxpos))];
    fdp_can_pos(rlz) = sum((adjcan (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(adjcan (idxpoc) <= q);
    fdp_com_pos(rlz) = sum((adjcom (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(adjcom (idxpoc) <= q);
    fdp_two_pos(rlz) = sum((adjtwo (idxpos) <= q) & ~ maskPos(idxpos)) / sum(adjtwo (idxpos) <= q);
    fdp_spl_pos(rlz) = sum((adjspl (idxpos) <= q) & ~ maskPos(idxpos)) / sum(adjspl (idxpos) <= q);
    fdp_Old_pos(rlz) = sum((adjOld (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(adjOld (idxpoc) <= q);
    fdp_Q11_pos(rlz) = sum((adjQ11 (idxpos) <= q) & ~ maskPos(idxpos)) / sum(adjQ11 (idxpos) <= q);
    fdp_Q12_pos(rlz) = sum((adjQ12 (idxpoc) <= q) & ~ maskPos(idxpos)) / sum(adjQ12 (idxpoc) <= q);

    % Empirical FDRs (negative side of the map, i.e., user interested in negative results only, i.e., results opposite to the direction of the contrast)
    idxnec = [false(size(idxneg));idxneg];
    fdp_can_neg(rlz) = sum((adjcan (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(adjcan (idxnec) <= q);
    fdp_com_neg(rlz) = sum((adjcom (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(adjcom (idxnec) <= q);
    fdp_two_neg(rlz) = sum((adjtwo (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(adjtwo (idxneg) <= q);
    fdp_spl_neg(rlz) = sum((adjspl (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(adjspl (idxneg) <= q);
    fdp_Old_neg(rlz) = sum((adjOld (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(adjOld (idxnec) <= q);
    fdp_Q11_neg(rlz) = sum((adjQ11 (idxneg) <= q) & ~ maskNeg(idxneg)) / sum(adjQ11 (idxneg) <= q);
    fdp_Q12_neg(rlz) = sum((adjQ12 (idxnec) <= q) & ~ maskNeg(idxneg)) / sum(adjQ12 (idxnec) <= q);
end

% If the FDP is NaN, then there were no positives, which we'd interpret as
% zero false positives (as opposed to 0/0 false positives)
fdp_can     (isnan(fdp_can)) = 0;
fdp_com     (isnan(fdp_com)) = 0;
fdp_two     (isnan(fdp_two)) = 0;
fdp_spl     (isnan(fdp_spl)) = 0;
fdp_Old     (isnan(fdp_Old)) = 0;
fdp_Q11     (isnan(fdp_Q11)) = 0;
fdp_Q12     (isnan(fdp_Q12)) = 0;
fdp_can_pos (isnan(fdp_can_pos)) = 0;
fdp_com_pos (isnan(fdp_com_pos)) = 0;
fdp_two_pos (isnan(fdp_two_pos)) = 0;
fdp_spl_pos (isnan(fdp_spl_pos)) = 0;
fdp_Old_pos (isnan(fdp_Old_pos)) = 0;
fdp_Q11_pos (isnan(fdp_Q11_pos)) = 0;
fdp_Q12_pos (isnan(fdp_Q12_pos)) = 0;
fdp_can_neg (isnan(fdp_can_neg)) = 0;
fdp_com_neg (isnan(fdp_com_neg)) = 0;
fdp_two_neg (isnan(fdp_two_neg)) = 0;
fdp_spl_neg (isnan(fdp_spl_neg)) = 0;
fdp_Old_neg (isnan(fdp_Old_neg)) = 0;
fdp_Q11_neg (isnan(fdp_Q11_neg)) = 0;
fdp_Q12_neg (isnan(fdp_Q12_neg)) = 0;

% Compute results
J.BothSides.canonical  = [mean(fdp_can),     confint(fdp_can)];
J.BothSides.combined   = [mean(fdp_com),     confint(fdp_com)];
J.BothSides.twotailed  = [mean(fdp_two),     confint(fdp_two)];
J.BothSides.split2tail = [mean(fdp_spl),     confint(fdp_spl)];
J.BothSides.old        = [mean(fdp_Old),     confint(fdp_Old)];
J.BothSides.question11 = [mean(fdp_Q11),     confint(fdp_Q11)];
J.BothSides.question12 = [mean(fdp_Q12),     confint(fdp_Q12)];
J.PosSide.canonical    = [mean(fdp_can_pos), confint(fdp_can_pos)];
J.PosSide.combined     = [mean(fdp_com_pos), confint(fdp_com_pos)];
J.PosSide.twotailed    = [mean(fdp_two_pos), confint(fdp_two_pos)];
J.PosSide.split2tail   = [mean(fdp_spl_pos), confint(fdp_spl_pos)];
J.PosSide.old          = [mean(fdp_Old_pos), confint(fdp_Old_pos)];
J.PosSide.question11   = [mean(fdp_Q11_pos), confint(fdp_Q11_pos)];
J.PosSide.question12   = [mean(fdp_Q12_pos), confint(fdp_Q12_pos)];
J.NegSide.canonical    = [mean(fdp_can_neg), confint(fdp_can_neg)];
J.NegSide.combined     = [mean(fdp_com_neg), confint(fdp_com_neg)];
J.NegSide.twotailed    = [mean(fdp_two_neg), confint(fdp_two_neg)];
J.NegSide.split2tail   = [mean(fdp_spl_neg), confint(fdp_spl_neg)];
J.NegSide.old          = [mean(fdp_Old_neg), confint(fdp_Old_neg)];
J.NegSide.question11   = [mean(fdp_Q11_neg), confint(fdp_Q11_neg)];
J.NegSide.question12   = [mean(fdp_Q12_neg), confint(fdp_Q12_neg)];
J.elapsed              = toc(st);

% Write results
writejson(J,outputfile)

function run_scenario(commonfile,scenariofile,outputfile)

% Keep track of time
st = tic;

% Read general configuration (common to all scenarios)
Jfdp = readjson(commonfile);

% Read json with scenario configuration
Js = readjson(scenariofile);

% Merge them, for the outputs later
F = fieldnames(Js);
for f = 1:numel(F)
    Jfdp.(F{f}) = Js.(F{f});
end

% Simulation parameters
if isoctave()
    rand('seed',Jfdp.seed) %#ok<RAND>
else
    rng(Jfdp.seed);
end
numTests        = Jfdp.numTests;        % number of "voxels"
effectSize      = Jfdp.effectSize;      % z-stat that determines the synthetic effect size
fracPos         = Jfdp.fracPos;         % fraction of tests with true positive effect
fracNeg         = Jfdp.fracNeg;         % fraction of tests with true negative effect
rho             = Jfdp.rho;             % sign of the compound symmetric correlation among the numTests
q               = Jfdp.q;               % test level, E(FDR) to be controlled
numRealizations = Jfdp.numRealizations; % number of times we repeat the simulation
FDRmethod       = Jfdp.FDRmethod;       % use 'bh1995' or 'bky2006'
CImethod        = Jfdp.CImethod;        % use 'Wald' or 'Wilson'
alpha           = Jfdp.alpha;           % for the confidence interval
Jpwr            = Jfdp;                 % to save power

% Vars for later
fdp_can     = zeros(numRealizations,1);
fdp_com     = zeros(numRealizations,1);
fdp_two     = zeros(numRealizations,1);
fdp_spl     = zeros(numRealizations,1);
fdp_can_pos = zeros(numRealizations,1);
fdp_com_pos = zeros(numRealizations,1);
fdp_two_pos = zeros(numRealizations,1);
fdp_spl_pos = zeros(numRealizations,1);
fdp_can_neg = zeros(numRealizations,1);
fdp_com_neg = zeros(numRealizations,1);
fdp_two_neg = zeros(numRealizations,1);
fdp_spl_neg = zeros(numRealizations,1);
pwr_can     = zeros(numRealizations,1);
pwr_com     = zeros(numRealizations,1);
pwr_two     = zeros(numRealizations,1);
pwr_spl     = zeros(numRealizations,1);
pwr_can_pos = zeros(numRealizations,1);
pwr_com_pos = zeros(numRealizations,1);
pwr_two_pos = zeros(numRealizations,1);
pwr_spl_pos = zeros(numRealizations,1);
pwr_can_neg = zeros(numRealizations,1);
pwr_com_neg = zeros(numRealizations,1);
pwr_two_neg = zeros(numRealizations,1);
pwr_spl_neg = zeros(numRealizations,1);

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
truePos = (1:numTests)'    <= numPos;
trueNeg = (numTests:-1:1)' <= numNeg;
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
    testPos = zstats > 0;
    testNeg = ~testPos;

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
    [~,~,adjspl(testPos)] = fdrfun(pvals2(testPos));
    [~,~,adjspl(testNeg)] = fdrfun(pvals2(testNeg));

    % Empirical FDRs (global, i.e., the user looks into both sides of the map)
    fdp_can(rlz) = sum((adjcan <= q) & ~ [truePos;trueNeg]) / sum(adjcan <= q);
    fdp_com(rlz) = sum((adjcom <= q) & ~ [truePos;trueNeg]) / sum(adjcom <= q);
    fdp_two(rlz) = sum((adjtwo <= q) & ~ (truePos|trueNeg)) / sum(adjtwo <= q);
    fdp_spl(rlz) = sum((adjspl <= q) & ~ (truePos|trueNeg)) / sum(adjspl <= q);

    % Empirical FDRs (positive side of the map, i.e., user interested in positive results only, i.e., results that match the direction of the contrast)
    testPoc = [testPos;false(size(testPos))];
    fdp_can_pos(rlz) = sum((adjcan (testPoc) <= q) & ~ truePos(testPos)) / sum(adjcan (testPoc) <= q);
    fdp_com_pos(rlz) = sum((adjcom (testPoc) <= q) & ~ truePos(testPos)) / sum(adjcom (testPoc) <= q);
    fdp_two_pos(rlz) = sum((adjtwo (testPos) <= q) & ~ truePos(testPos)) / sum(adjtwo (testPos) <= q);
    fdp_spl_pos(rlz) = sum((adjspl (testPos) <= q) & ~ truePos(testPos)) / sum(adjspl (testPos) <= q);

    % Empirical FDRs (negative side of the map, i.e., user interested in negative results only, i.e., results opposite to the direction of the contrast)
    testNec = [false(size(testNeg));testNeg];
    fdp_can_neg(rlz) = sum((adjcan (testNec) <= q) & ~ trueNeg(testNeg)) / sum(adjcan (testNec) <= q);
    fdp_com_neg(rlz) = sum((adjcom (testNec) <= q) & ~ trueNeg(testNeg)) / sum(adjcom (testNec) <= q);
    fdp_two_neg(rlz) = sum((adjtwo (testNeg) <= q) & ~ trueNeg(testNeg)) / sum(adjtwo (testNeg) <= q);
    fdp_spl_neg(rlz) = sum((adjspl (testNeg) <= q) & ~ trueNeg(testNeg)) / sum(adjspl (testNeg) <= q);

    % Empirical power (global, i.e., user looks into both sides of the map)
    pwr_can(rlz) = sum(adjcan <= q) / sum([truePos;trueNeg]);
    pwr_com(rlz) = sum(adjcom <= q) / sum([truePos;trueNeg]);
    pwr_two(rlz) = sum(adjtwo <= q) / sum((truePos|trueNeg));
    pwr_spl(rlz) = sum(adjspl <= q) / sum((truePos|trueNeg));

    % Empirical power (positive side of the map)
    pwr_can_pos(rlz) = sum(adjcan (testPoc) <= q) / sum(truePos(testPos));
    pwr_com_pos(rlz) = sum(adjcom (testPoc) <= q) / sum(truePos(testPos));
    pwr_two_pos(rlz) = sum(adjtwo (testPos) <= q) / sum(truePos(testPos));
    pwr_spl_pos(rlz) = sum(adjspl (testPos) <= q) / sum(truePos(testPos));

    % Empirical power (negative side of the map)
    pwr_can_neg(rlz) = sum(adjcan (testNec) <= q) / sum(trueNeg(testNeg));
    pwr_com_neg(rlz) = sum(adjcom (testNec) <= q) / sum(trueNeg(testNeg));
    pwr_two_neg(rlz) = sum(adjtwo (testNeg) <= q) / sum(trueNeg(testNeg));
    pwr_spl_neg(rlz) = sum(adjspl (testNeg) <= q) / sum(trueNeg(testNeg));
end

% If the FDP is NaN, then there were no positives, which we'd interpret as
% zero false positives (as opposed to 0/0 false positives)
fdp_can     (isnan(fdp_can)) = 0;
fdp_com     (isnan(fdp_com)) = 0;
fdp_two     (isnan(fdp_two)) = 0;
fdp_spl     (isnan(fdp_spl)) = 0;
fdp_can_pos (isnan(fdp_can_pos)) = 0;
fdp_com_pos (isnan(fdp_com_pos)) = 0;
fdp_two_pos (isnan(fdp_two_pos)) = 0;
fdp_spl_pos (isnan(fdp_spl_pos)) = 0;
fdp_can_neg (isnan(fdp_can_neg)) = 0;
fdp_com_neg (isnan(fdp_com_neg)) = 0;
fdp_two_neg (isnan(fdp_two_neg)) = 0;
fdp_spl_neg (isnan(fdp_spl_neg)) = 0;

% This is unnecessary for power but let's do anyway for symmetry
pwr_can     (isnan(pwr_can)) = 0;
pwr_com     (isnan(pwr_com)) = 0;
pwr_two     (isnan(pwr_two)) = 0;
pwr_spl     (isnan(pwr_spl)) = 0;
pwr_can_pos (isnan(pwr_can_pos)) = 0;
pwr_com_pos (isnan(pwr_com_pos)) = 0;
pwr_two_pos (isnan(pwr_two_pos)) = 0;
pwr_spl_pos (isnan(pwr_spl_pos)) = 0;
pwr_can_neg (isnan(pwr_can_neg)) = 0;
pwr_com_neg (isnan(pwr_com_neg)) = 0;
pwr_two_neg (isnan(pwr_two_neg)) = 0;
pwr_spl_neg (isnan(pwr_spl_neg)) = 0;

% Compute results
Jfdp.BothSides.canonical  = [mean(fdp_can),     confint(fdp_can)];
Jfdp.BothSides.combined   = [mean(fdp_com),     confint(fdp_com)];
Jfdp.BothSides.twotailed  = [mean(fdp_two),     confint(fdp_two)];
Jfdp.BothSides.split2tail = [mean(fdp_spl),     confint(fdp_spl)];
Jfdp.PosSide.canonical    = [mean(fdp_can_pos), confint(fdp_can_pos)];
Jfdp.PosSide.combined     = [mean(fdp_com_pos), confint(fdp_com_pos)];
Jfdp.PosSide.twotailed    = [mean(fdp_two_pos), confint(fdp_two_pos)];
Jfdp.PosSide.split2tail   = [mean(fdp_spl_pos), confint(fdp_spl_pos)];
Jfdp.NegSide.canonical    = [mean(fdp_can_neg), confint(fdp_can_neg)];
Jfdp.NegSide.combined     = [mean(fdp_com_neg), confint(fdp_com_neg)];
Jfdp.NegSide.twotailed    = [mean(fdp_two_neg), confint(fdp_two_neg)];
Jfdp.NegSide.split2tail   = [mean(fdp_spl_neg), confint(fdp_spl_neg)];
Jfdp.elapsed              = toc(st);
Jpwr.BothSides.canonical  = [mean(pwr_can),     confint(pwr_can)];
Jpwr.BothSides.combined   = [mean(pwr_com),     confint(pwr_com)];
Jpwr.BothSides.twotailed  = [mean(pwr_two),     confint(pwr_two)];
Jpwr.BothSides.split2tail = [mean(pwr_spl),     confint(pwr_spl)];
Jpwr.PosSide.canonical    = [mean(pwr_can_pos), confint(pwr_can_pos)];
Jpwr.PosSide.combined     = [mean(pwr_com_pos), confint(pwr_com_pos)];
Jpwr.PosSide.twotailed    = [mean(pwr_two_pos), confint(pwr_two_pos)];
Jpwr.PosSide.split2tail   = [mean(pwr_spl_pos), confint(pwr_spl_pos)];
Jpwr.NegSide.canonical    = [mean(pwr_can_neg), confint(pwr_can_neg)];
Jpwr.NegSide.combined     = [mean(pwr_com_neg), confint(pwr_com_neg)];
Jpwr.NegSide.twotailed    = [mean(pwr_two_neg), confint(pwr_two_neg)];
Jpwr.NegSide.split2tail   = [mean(pwr_spl_neg), confint(pwr_spl_neg)];
Jpwr.elapsed              = toc(st);

% Write results
writejson(Jfdp,strrep(outputfile,'.json','_fdp.json'));
writejson(Jpwr,strrep(outputfile,'.json','_pwr.json'));


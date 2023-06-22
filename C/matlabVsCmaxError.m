%compare performance of C and Matlab bky implementations
% Simulation parameters
rng('shuffle');
rng(3)
numTests   = 10000;
effectSize = 3;
fracPos    = 0.2;
fracNeg    = 0.4;
rho        = 0.8;  % Compound symmetric correlation among the numTests (rho>=0)
q          = 0.05; % test level, E(FDR) to be controlled
nRlz       = 1;
alpha      = 0.05; % for the confidence interval
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
maxDifference = 0.0;
maxDifferenceAdj = 0.0;
for rlz = 1:nRlz

    % Create random data, add signal
    zstats = signal + sqrt(1-rho)*randn(numTests,1) + sqrt(rho)*randn(1,1);
    idxpos = zstats > 0;
    idxneg = ~idxpos;

    % The usual pvalues, one-tailed
    pvals = normcdf(zstats,'upper');
    pvals = double(pvals);
    %test with Matlab
    [pthr,~,padj] = bky7(pvals);

    %compare with C code:
    pth = fileparts(mfilename('fullpath'));
    exefnm = fullfile(pth, 'bky');
    pvalfnm = fullfile(pth, 'pvals.bin');
    padjfnm = fullfile(pth, 'padj.bin');

    %write binary file to disk
    fileID = fopen(pvalfnm,'wb');
    fwrite(fileID, double(pvals),'double');
    fclose(fileID);
    
    %test with C executable
    if (~exist(exefnm,"file"))
        error("Unable to find executable %s", exefnm)
    end
    cmd = sprintf('%s %s 0.05 n %s', exefnm, pvalfnm, padjfnm);
    [status,cmdout]=system(cmd);
    cthr = str2double(cmdout);
    %read C padj
    fileID = fopen(padjfnm);
    cpadj = fread(fileID,[numTests 1],'double');
    fclose(fileID);

    %compare Matlab vs C
    diff = abs(pthr - cthr);
    maxDifference = max(maxDifference, diff);
    diffAdj = max(abs(padj(:) - cpadj(:)));
    maxDifferenceAdj = max(maxDifferenceAdj, diffAdj);
    
    dx = abs(padj(:) - cpadj(:));
    [mx, idx] = max(dx)
    padj(idx)
    cpadj(idx)
end
if (maxDifference == 0.0)
    fprintf("Identical results\n");
else
    fprintf("Maximum difference %g\n", maxDifference);
end
if (maxDifferenceAdj == 0.0)
    fprintf("Identical pAdj results\n");
else
    fprintf("Maximum pAdj difference %g\n", maxDifferenceAdj);
end
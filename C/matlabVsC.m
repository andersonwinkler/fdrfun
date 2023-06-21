%compare performance of C and Matlab bky implementations
% Simulation parameters
rng('shuffle');
numTests   = 10000;
effectSize = 3;
fracPos    = 0.2;
fracNeg    = 0.4;
rho        = 0.8;  % Compound symmetric correlation among the numTests (rho>=0)
q          = 0.05; % test level, E(FDR) to be controlled
nRlz       = 10000;
alpha      = 0.05; % for the confidence interval
% For each realization
maxDifference = 0.0;
for rlz = 1:nRlz

    % Create random data, add signal
    zstats = signal + sqrt(1-rho)*randn(numTests,1) + sqrt(rho)*randn(1,1);
    idxpos = zstats > 0;
    idxneg = ~idxpos;

    % The usual pvalues, one-tailed
    pvals = normcdf(zstats,'upper');
    pvals = single(pvals);
    %test with Matlab
    pthr = bky7(pvals);

    %write binary file to disk
    fileID = fopen('singles.bin','w');
    fwrite(fileID, single(pvals),'single');
    fclose(fileID);

    %test with C executable
    pth = fileparts(mfilename('fullpath'));
    exe = fullfile(pth, 'bky');
    if (~exist(exe,"file"))
        error("Unable to find executable %s", exe)
    end
    [status,cmdout]=system(sprintf('%s singles.bin 0.05 n', exe));
    cthr = str2double(cmdout);

    diff = abs(pthr - cthr);
    maxDifference = max(maxDifference, diff);

end
if (maxDifference == 0.0)
    fprintf("Identical results\n");
else
    fprintf("Maximum difference %g\n", maxDifference);
end
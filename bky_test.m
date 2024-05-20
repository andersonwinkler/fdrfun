% Raw p-values, already sorted in ascending order
p = [0.0026, 0.01, 0.014, 0.025, 0.042, 0.066, 0.1, 0.12, 0.17, 0.28, 0.36, 0.524, 0.61, 0.68, 0.78, 0.9, 0.96];

% Test level. The q = 0.08 is one in which we see a discrepancy between
% multiple.down and what the paper would have given, as demonstrated by the
% inefficient ("brute force") code below
q = 0.08;

% Number of tests
m = numel(p);

% The initial k; this will be replaced later, or not if no tests are found
% significant
k = 0;

% Also, let's start by declaring all hypotheses false
signif = false(1,m);

% It's a step-up procedure. Since we we look for the
% max i that satisfies a condition, we start from m and move downwards,
% not from 1 and upwards:
for i = m:-1:1

    % The "exists" here has the same role as the "exists" in step 1 of
    % Definition 7 of the BKY paper (page 496). We'll populate this
    % variable further down
    exists = false(1,i);
    
    % We'll check that that condition is satisfied for all j <= i
    for j = 1:i

        % We'll check that there exists at least one l >= j in which
        % that condition is satisfied 
        for l = j:m

            % If the condition is satisfied for at least one l, there is no
            % need to continue, so we mark that as satisfied and break the
            % loop so we can go to the next j
            if p(l) <= q*l/(m+1-j*(1-q))
                exists(j) = true;
                break
            end
        end
    end

    % Since we started from m, moving dowards towards the first p-values, then
    % then if for the current i, the condition is satisfied for all j <= 1
    % then we can stop here and call this i our critical k, since this is
    % the maximum i. Otherwise, we continue with the next smaller i.
    if all(exists)
        k = i;
        break
    end
end

% Now we declare the tests 1:k as significant. If no such k was found, then
% k is still 0, and all hypotheses remain declared false
signif(1:k) = true;

[repodir,~,~] = fileparts(mfilename('fullpath'));

% Load the table, sort
T = readtable(fullfile(repodir,'summary_results.csv'));
T = sortrows(T,{'FDRmethod','scenario'},{'ascend','ascend'});

% Select columns of interest, drop negative dep scenarios, convert to percentages
nidx = T.scenario <= 10;
Tx = T(nidx,14:end-1);
Tx = table2array(Tx)*100;

% Column headers
Th = T.Properties.VariableNames;
Th = Th(14:end-1);
hidx = contains(Th,'_question') | contains(Th,'_old');
Th = Th(:,hidx);

% Rearrange the table
Tx = Tx(:,hidx);
Ty = [];
Ty( 1:10, 1:9 ) = Tx( 1:10, 1:9);
Ty( 1:10,10:18) = Tx(11:20, 1:9);
Ty(11:20, 1:9 ) = Tx( 1:10,10:18);
Ty(11:20,10:18) = Tx(11:20,10:18);
Ty(21:30, 1:9 ) = Tx( 1:10,19:27);
Ty(21:30,10:18) = Tx(11:20,19:27);
for r = 1:size(Ty,1)
    if r > 20
        scn = r-20;
    elseif r > 10
        scn = r-10;
    else
        scn = r;
    end
    fprintf('\\textsc{%s}',lower(integer2roman(scn)));
    fprintf(' & %0.1f \\scalebox{.7}[1.0]{(%0.1f--%0.1f)}', Ty(r,:));
    fprintf(' \\\\\n');
end

% =========================================================================
function roman = integer2roman(n)
if (n < 1 || n > 3999)
    error('Input must be an integer between 1 and 3999.');
end
val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
syms = {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
roman = '';
for i = 1:length(val)
    while n >= val(i)
        n = n - val(i);
        roman = [roman, syms{i}];
    end
end
end

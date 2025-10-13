[repodir,~,~] = fileparts(mfilename('fullpath'));
noCI = false;
% Load the table, sort
T = readtable(fullfile(repodir,'summary_results_pwr.csv'));
T = sortrows(T,{'FDRmethod','scenario'},{'ascend','ascend'});

% Select columns of interest, convert to percentages
Tx = T(:,14:end-1);
Tx = table2array(Tx)*100;

% % Way 1: BH at the top, BKY at the bottom
% for r = 1:size(Tx,1)
%     if r > 10
%         scn = r-10;
%     else
%         scn = r;
%     end
%     fprintf('\\textsc{%s}',lower(integer2roman(scn)));
%     fprintf(' & %0.1f \\scalebox{.7}[1.0]{(%0.1f--%0.1f)}', Tx(r,:));
%     fprintf(' \\\\\n');
% end

% Way 2: BH on the left, BKY on the right
% Rearrange the table
Ty = [];
Ty( 1:10, 1:18) = Tx( 1:10,1:18);
Ty( 1:10,19:36) = Tx(11:20,1:18);
Ty(11:20, 1:18) = Tx( 1:10,19:36);
Ty(11:20,19:36) = Tx(11:20,19:36);
Ty(21:30, 1:18) = Tx( 1:10,37:54);
Ty(21:30,19:36) = Tx(11:20,37:54);

if noCI
    Ty = Ty(:,1:3:end);
end

for r = 1:size(Ty,1)
    if r > 20
        scn = r-20;
    elseif r > 10
        scn = r-10;
    else
        scn = r;
    end

    if noCI
        fprintf('%s',upper(integer2roman(scn)));
        fprintf(' & %0.1f', Ty(r,:));
    else
        fprintf('\\textsc{%s}',lower(integer2roman(scn)));
        fprintf(' & %0.1f \\scalebox{.7}[1.0]{(%0.1f--%0.1f)}', Ty(r,:));
    end
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

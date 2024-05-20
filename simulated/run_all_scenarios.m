% Directory of this repo
[repodir,~,~] = fileparts(mfilename('fullpath'));
addpath(repodir);
addpath(fullfile(repodir,'..'));

% List of configuration files defining the scenarios
Jfiles = dir(fullfile(repodir,'configs','*.json'));

% Run each of these scenarios
for j = 1:numel(Jfiles)
    fprintf('Running scenario %d/%d\n',j,numel(Jfiles));
    run_scenario_reviewer(...
        fullfile(repodir,'configs',Jfiles(j).name),...
        fullfile(repodir,'common_config.json'),...
        fullfile(repodir,'results',Jfiles(j).name));
end

% Generate a table with the results, save it
T = table();
warning('off','MATLAB:table:RowsAddedExistingVars');
for j = 1:numel(Jfiles)
    resultsfile = fullfile(repodir,'results',Jfiles(j).name);
    [~,scenario,~] = fileparts(resultsfile);
    J = readjson(resultsfile);
    F = fieldnames(J);
    for f = 1:numel(F)
        if isstruct(J.(F{f}))
            Fs = fieldnames(J.(F{f}));
            for fs = 1:numel(Fs)
                T{scenario,sprintf('%s_%s_mean', F{f},Fs{fs})} = J.(F{f}).(Fs{fs})(1);
                T{scenario,sprintf('%s_%s_lower',F{f},Fs{fs})} = J.(F{f}).(Fs{fs})(2);
                T{scenario,sprintf('%s_%s_upper',F{f},Fs{fs})} = J.(F{f}).(Fs{fs})(3);
            end
        elseif ischar(J.(F{f}))
            T{scenario,F{f}} = {J.(F{f})};
        else
            T{scenario,F{f}} = J.(F{f});
        end
    end
end
warning('on','MATLAB:table:RowsAddedExistingVars');
writetable(T,fullfile(repodir,'summary_results.csv'),'WriteRowNames',true);
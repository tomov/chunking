% model_exp_v2_3 but for ISL
%
clear all;
close all;

h = init_hyperparams;

start = [6 6 6 6 6 6];
goal = [1 1 1 1 1 1];
%ordinal = [1 2 3 4 5]; <-- don't use that, e.g. if we happened to have 6 1 by chance
index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..; ORDER CRUCIAL
% IMPORTANT first must be HBFS preferred action
nexts = [
7 5;
7 5;
7 5;
7 5;
7 5;
7 5
];

H = [1 1 1 2 2 2 3 3 3 3]; % ground truth


prefix = 'pf4'; % isl_MH or pf2
cachefile = 'isl_fig2_pf4.mat';
%prefix = 'isl_MH'; % isl_MH or pf2

templates = {sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=34.mat', prefix), 
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=68.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=103.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=150.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=197.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=246.mat', prefix)};
filename = sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', prefix);

[~, name] = system('hostname');
if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
    % local
    templates
    filename

else
    % cannon
    for t = 1:length(templates)
        templates{t} = fullfile(getenv('MY_SCRATCH'), 'chunking', templates{t});
    end
    filename = fullfile(getenv('MY_SCRATCH'), 'chunking', filename);

    templates
    filename

end


for t = 1:length(templates)

    files = dir(templates{t});

    load(filename, 'D');

    s = start(t);
    g = goal(t);

    for i = 1:length(files)
        files(i).name

        load(fullfile(files(i).folder, files(i).name), 'subj', 'particles');

        assert(s == D(subj).tasks.s(index(t)));
        assert(g == D(subj).tasks.g(index(t)));

        j = randi(length(particles)); % alternative: mode, or average

        [path, hpath] = hbfs(s, g, particles(j).H, D(subj));
        move(subj, t) = path(2);

        % eps-greedy: choose random neighbor w/ small prob 
        if rand() < 1 - h.eps
            move(subj, t) = datasample(find(D(subj).G.E(s,:)), 1);
        end
    end
    
end


cachefile
save(cachefile);

% load(cachefile);


mv = move == 5; %nexts(i,1);

ms = mean(mv, 1);
sems = std(mv, 1) / sqrt(size(mv, 1));

% swap to be consistent with other plots

filename
save(filename);

figure;
hold on;
bar(ms);
errorbar(ms, sems, 'linestyle', 'none', 'color', 'black');
plot([0 6], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 0.7], '-', 'color', [0.5 0.5 0.5])
hold off;
ylabel('p(HBFS direction)');
xticks(1:6);
xticklabels(index);
xlabel('trial #');
title(sprintf('model N = %d', size(move,1)));


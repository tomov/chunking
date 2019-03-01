
clear all;
close all;
%[data, Ts] = load_data('exp/results', 246); % for exp_v2_3 (subway 10 unlearn)

h = init_hyperparams;
nsamples = 100;
burnin = 1;
lag = 1;

% from analyze_data
% for exp_v2_3 subway 10 unlearn 
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

% from model_all_data
D = init_Ds_from_data('exp/results/exp_v2_3_subway10_unlearn', true);
D_full = D;
% TODO rm & uncomment above
save tmp.mat;
%load tmp.mat;

for subj = 1:length(D) % for each subject
    D(subj).tasks.s = [];
    D(subj).tasks.g = [];
    for i = 1:D(subj).G.N
        D(subj).r{i} = [];
    end
    subj

    [samples, post] = sample(D(subj), h, nsamples, burnin, lag);
    H(subj) = samples(end);

    t = 1;
    for i = 1:length(index)
        i
        while t < index(i)
            D(subj).tasks.s = [D(subj).tasks.s D_full(subj).tasks.s(t)];
            D(subj).tasks.g = [D(subj).tasks.g D_full(subj).tasks.g(t)];
            t = t + 1;
        end
        
        s = start(i);
        g = goal(i);
        assert(s == D_full(subj).tasks.s(t));
        assert(g == D_full(subj).tasks.g(t));
        
        [samples, post] = sample(D(subj), h, nsamples, burnin, lag, H(subj));
        H(subj) = samples(end);

        [path, hpath] = hbfs(s, g, H(subj), D(subj));
        move(subj, i) = path(2) == nexts(i,1);
    end
end

save('model_exp_v2_3.mat');

ms = mean(move, 1);
sems = std(move, 1) / sqrt(size(move, 1));

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
title(sprintf('model N = %d', length(data)));

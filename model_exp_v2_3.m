
clear all;
close all;


h = init_hyperparams;
nsamples = 10000;
take_map = false;

sem = @(x) std(x) / sqrt(length(x));

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

assert(take_map == false); % don't support MAP
filename = sprintf('model_exp_v2_3_circ_alpha=%.4f_nsamples=%d_div_eps=%.4f_last.mat', h.alpha, nsamples, h.eps);
filename

%{
D = init_Ds_from_data('exp/results/exp_v2_3_subway10_unlearn_circ', true);
D_full = D;

tic

for subj = 1:length(D) % for each subject
    D(subj).tasks.s = [];
    D(subj).tasks.g = [];
    for i = 1:D(subj).G.N
        D(subj).r{i} = [];
    end

    fprintf('infer H: subject %d\n', subj);

    H(subj) = sample_c(D(subj), h, 1); % take 1 sample to start with

    t = 1;
    for i = 1:length(index)
        while t < index(i)
            D(subj).tasks.s = [D(subj).tasks.s D_full(subj).tasks.s(t)];
            D(subj).tasks.g = [D(subj).tasks.g D_full(subj).tasks.g(t)];
            t = t + 1;
        end
        
        H(subj) = sample_c(D(subj), h, 1, round(nsamples / length(index)), 1, H(subj));

        chosen_H{subj, i} = H(subj);
    end
end

toc

save(filename, '-v7.3');
%}
load(filename);

clear move;
for subj = 1:length(D) % for each subject
    for i = 1:length(index)

        fprintf('HBFS: subject %d\n', subj);

        s = start(i);
        g = goal(i);
        assert(s == D_full(subj).tasks.s(index(i)));
        assert(g == D_full(subj).tasks.g(index(i)));
       
        [path, hpath] = hbfs(s, g, chosen_H{subj, i}, D(subj));
        move(subj, i) = path(2);

        % eps-greedy: choose random neighbor w/ small prob 
        if rand() < 1 - h.eps
            move(subj, i) = datasample(find(D(subj).G.E(s,:)), 1);
        end
    end
end

mv = move == nexts(i,1);

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



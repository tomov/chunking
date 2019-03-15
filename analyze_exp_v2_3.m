% analyze behavioral data from subway 10 unlearn with circles

clear all;
close all;

[data, Ts] = load_data('exp/results/exp_v2_3_subway10_unlearn_circ', 246, false); % for exp_v2_3 (subway 10 unlearn)


sem = @(x) std(x) / sqrt(length(x));

s = [];
g = [];
len = [];
group = [];
dir = []; % direction = 2nd state on path
ord = []; % ordinal of trial type within phase (e.g. "first 1->6", "second 1->6", etc)
subj_group = [];
subj_len = [];
s_id = []; % subject index
t_id = []; % trial index
for subj = 1:size(data,1) % for each subject
    %phase = 1; % = 1 for training exp_v3_7, usually it's 2 = test; TODO DON'T FORGET TO CHANGE TO 2 for subway10 and other old stuff
    for phase = 1:2 % for exp_v2_3!!! TODO undo
        for i = 1:length(data(subj, phase).s) % for each trial 
            which = find(data(subj, phase).s == data(subj, phase).s(i) & data(subj, phase).g == data(subj, phase).g(i));
            clear o;
            o(which) = find(which);
            ord = [ord; o(i)];
            s = [s; data(subj, phase).s(i)];
            g = [g; data(subj, phase).g(i)];
            len = [len; data(subj, phase).len(i)];
            dir = [dir; data(subj, phase).path{i}(2)];
            group = [group; data(subj, phase).group(i)];
            s_id = [s_id; subj];
            t_id = [t_id; i + (phase - 1) * 103]; % TODO hack for exp_v2_3
        end
    end
    subj_group = [subj_group; data(subj,1).group(1)];
    subj_len = [subj_len; mean(data(subj, 1).len)];
end


% show learning
%

figure;
ms = [];
es = [];
for t = 1:length(data(1,1).len)
    l = [];
    for subj = 1:size(data,1)
        l = [l data(subj,1).len(t)];
    end
    m = mean(l);
    e = std(l) / sqrt(length(l));
    ms = [ms m];
    es = [es e];
end

subplot(2,1,1);
errorbar(ms, es);
xlabel('training trial');
ylabel('path length');
title('all trials');



% for exp_v2_3 subway 10 unlearn 
start = [6 6 6 6 6 6];
goal = [1 1 1 1 1 1];
%ordinal = [1 2 3 4 5]; <-- don't use that, e.g. if we happened to have 6 1 by chance
index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..
nexts = [
7 5;
7 5;
7 5;
7 5;
7 5;
7 5
];



figure;

pl.tests = [3 3 3 3 3 3]; % 1 = right tailed, 2 = left tailed, 3 = two-tailed


ms = [];
sems = [];
for t = 1:length(start)
    %which = s == start(t) & g == goal(t) & ord == ordinal(t); % ord is usually just 1
    which = s == start(t) & g == goal(t) & t_id == index(t); % ord is usually just 1
    move = dir(which);
    m = nexts(t,:);
    c1 = sum(move == m(1)); % count 1
    c2 = sum(move == m(2)); % count 2
    n = sum(which);

    switch pl.tests(t)
        case 1 % right-tailed
            p = 1 - binocdf(c1, n, 0.5);
        case 2 % left-tailed
            p = binocdf(c1, n, 0.5);
        case 3 % two-tailed
            p = 2 * binocdf(min(c1,c2), n, 0.5);
        otherwise
            assert(false);
    end

    ms = [ms mean(move == m(1))];
    sems = [sems sem(move == m(1))];

    subplot(2,3,t);
    bar(1:2, [c1 c2]);
    hold on;
    y = binoinv([0.025 0.975], n, 0.5);


    pl.ci(t) = (y(2) - y(1)) / 2;
    pl.n(t) = n;
    pl.m(t) = c1;
    pl.p(t) = p;

    plot([0 3], [y(1) y(1)], '--', 'Color', [0.5 0.5 0.5]);
    plot([0 3], [y(2) y(2)], '--', 'Color', [0.5 0.5 0.5]);
    hold off;
    xticklabels({num2str(m(1)), num2str(m(2))});
    title(sprintf('%d -> %d: p = %.3f (c1 = %d, n = %d)', start(t), goal(t), p, c1, n));

    fprintf('trial #%d (%d -> %d): p = %.3f (c1 = %d, n = %d)\n', index(t), start(t), goal(t), p, c1, n);
end

figure;
hold on;
bar(ms);
errorbar(ms, sems, 'linestyle', 'none', 'color', 'black');
plot([0 6], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 0.7], '-', 'color', [0.5 0.5 0.5])
hold off;
assert(nexts(1,1) == 7);
ylabel('p(go to 7)');
xticks(1:6);
xticklabels(index);
xlabel('trial # (probe)');
title(sprintf('human N = %d', length(data)));


% compare 3rd probe and 4th probe trial
%
which_g1 = t_id == index(3);
which_g2 = t_id == index(4);
assert(nexts(3,1) == nexts(4,1));
[h, p, ci, stats] = ttest2(dir(which_g1) == nexts(3,1), dir(which_g2) == nexts(4,1));
p
stats


% stats for LEARNING -- is there a ramp?
%
assert(nexts(3,1) == 7);
assert(nexts(4,1) == 7);
direction = dir == 7;
[~,trial_idx] = max(t_id == index, [], 2); % convert to 1..3 
subject = s_id;
% subset
which_probes = ismember(t_id, index(1:3));
direction = direction(which_probes);
trial_idx = trial_idx(which_probes);
subject = subject(which_probes);
% create table
tbl = table(direction, trial_idx, subject);

formula = 'direction ~ 1 + trial_idx + (1 + trial_idx | subject)';
result1 = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result1);

H = [0 1];
[p, F, DF1, DF2] = coefTest(result1, H);
fprintf('Learning: is there a ramp on probe trials 1..3? coef for trial_idx = %f, p = %f, F(%d,%d) = %f\n', H * beta, p, DF1, DF2, F);



% stats for UNLEARNing: trials 4..6
%
assert(nexts(3,1) == 7);
assert(nexts(4,1) == 7);
direction = dir == 7;
[~,trial_idx] = max(t_id == index, [], 2); % convert to 1..3 
subject = s_id;
% subset
which_probes = ismember(t_id, index(4:6));
direction = direction(which_probes);
trial_idx = trial_idx(which_probes);
subject = subject(which_probes);
% create table
tbl = table(direction, trial_idx, subject);

formula = 'direction ~ 1 + trial_idx + (1 + trial_idx | subject)';
result2 = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result2);

H = [1 0];
[p, F, DF1, DF2] = coefTest(result2, H);
fprintf('Unlearning: are probe trials 4..6 below 0.5? intercept = %f, p = %f, F(%d,%d) = %f\n', H * beta, p, DF1, DF2, F);


% stats for UNLEARNing: trials 3 and 4
%
assert(nexts(3,1) == 7);
assert(nexts(4,1) == 7);
direction = dir == 7;
[~,trial_idx] = max(t_id == index, [], 2); % convert to 1..3 
subject = s_id;
% subset
which_probes = ismember(t_id, index(3:4));
direction = direction(which_probes);
trial_idx = trial_idx(which_probes);
subject = subject(which_probes);
% create table
tbl = table(direction, trial_idx, subject);

formula = 'direction ~ 1 + trial_idx + (1 + trial_idx | subject)';
result3 = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result3);

H = [0 1];
[p, F, DF1, DF2] = coefTest(result3, H);
fprintf('Unlearning: is probe trial 4 below probe trial 3? coef for trial_idx = %f, p = %f, F(%d,%d) = %f\n', H * beta, p, DF1, DF2, F);

save('analyze_exp_v2_3.mat');

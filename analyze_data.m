% analyze behavioral data from chunking experiment

%[data, Ts] = load_data('exp/results', 165); % for exp_v3_7 (mail delivery map aka exp 3 scaled up)
%[data, Ts] = load_data('exp/results', 105); % for exp_v3_8 (subway 18 map aka mail delivery scaled down)
%[data, Ts] = load_data('exp/results', 81); % for exp_v1_6 (subway 10 but no assoc)
%[data, Ts] = load_data('exp/results', 101); % for exp_v2_1 (subway 10 no adj, no assoc)
%[data, Ts, ~, durs] = load_data('exp/results', 205); % for exp_v2_2 (subway 18 no adj, no assoc)
%[data, Ts, ~, durs] = load_data('exp/results/exp_v2_2_subway18_noadj_noassoc', 205, false); % for exp_v2_2 (subway 18 no adj, no assoc)
%[data, Ts] = load_data('exp/results/exp_v2_1_subway10_noadj_noassoc/', 101); % for exp_v2_1 (subway 10 no adj, no assoc)
%[data, Ts] = load_data('exp/results/ARCHIVE/exp_v2_1_batch_2/', 101); % for exp_v2_1 (subway 10 no adj, no assoc)
%[data, Ts] = load_data('exp/results/subway10_repro', 83); % for subway 10 repro TODO change phase = 2!
%[data, Ts] = load_data('exp/results/subway9', 81); % for subway 10 repro TODO change phase = 2!

%[data, Ts] = load_data('exp/results/exp_v2_3_subway10_unlearn/', 246); % for exp_v2_3 (subway 10 unlearn)
%[data, Ts] = load_data('exp/results/exp_v2_3_subway10_unlearn_circ', 246, false); % for exp_v2_3 (subway 10 unlearn)
%[data, Ts] = load_data('exp/results/ARCHIVE/exp_v2_3_subway10_unlearn_batch2', 246, false); % for exp_v2_3 (subway 10 unlearn)
%[data, Ts] = load_data('exp/results/exp_', 246, false); % for exp_v2_3 (subway 10 unlearn)
%[data, Ts] = load_data('exp/results/mines10_map', 101, false); % for exp_v2_3 (subway 10 unlearn)
%[data, Ts] = load_data('exp/results/exp_v2_1_subway10_noadj_noassoc', 101, false); 

load data.mat

%data = data(durs < 50, :);

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
    phase = 1; % = 1 for training exp_v3_7, usually it's 2 = test; TODO DON'T FORGET TO CHANGE TO 2 for subway10 and other old stuff
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

        t_id = [t_id; i];
        %t_id = [t_id; i + (phase - 1) * 103]; % TODO hack for exp_v2_3
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


% show learning on non-task trials

% TODO it's hardcoded for exp_v1_6
%len = [];
%for subj = 1:size(data,1)
%    which = ~ismember(data(subj,1).s, [2 4 10]) | ~ismember(data(subj,1).g, [3 5 7]);
%    l = data(subj,1).len(which);
%    l = [l nan(1, 41 - length(l))];
%    len = [len; l];
%end
%ms = mean(len, 1);
%sems = std(len, 1) / sqrt(size(len, 1));
%
%subplot(2,1,2);
%errorbar(ms, sems);
%xlabel('training trial');
%ylabel('path length');
%title('random trials');



% show test choices
%

% for exp_v2_1.html (subway 10 no adj no assoc)
start = [6];
goal = [1];
%ordinal = [1];
index = [101]; % from html -- @ ..
nexts = [
5 7
];


% for exp_v4.html (mines 10 map)
%{
start = [6];
goal = [1];
ordinal = [1];
nexts = [
5 7
];
%}

% for exp_v1_6.html (subway 10 no assoc), exp_v2_1.html
%start = [6];
%goal = [1];
%ordinal = [1];
%nexts = [
%7 5
%];

% for exp_v2_3 subway 10 unlearn 
%start = [6 6 6 6 6 6];
%goal = [1 1 1 1 1 1];
%%ordinal = [1 2 3 4 5]; <-- don't use that, e.g. if we happened to have 6 1 by chance
%index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..
%nexts = [
%7 5;
%7 5;
%7 5;
%7 5;
%7 5;
%7 5
%];


% for mail delivery  exp_v3_7.html
%start = [105 105 105 105 105];
%goal = [114 114 114 114 114];
%ordinal = [1 2 3 4 5];
%nexts = [
%106 104;
%106 104;
%106 104;
%106 104;
%106 104
%];

% for subway 10 and 9
%start = [6 7 3 1 2 8];
%goal = [1 2 8 6 7 3];
%nexts = [
%5 7;
%8 6;
%2 4;
%2 10;
%1 3;
%9 7];

% for subway 8
%{
start = [5 3 6 1 7 2];
goal = [1 7 2 5 3 6];
nexts = [
4 6;
2 4;
7 5;
2 8;
8 6;
1 3
];
%}

% for subway 6
%{
start = [4 5 2];
goal = [1 2 5];
nexts = [
3 5;
6 4;
1 3];
%}

% for subway 12
%
%{
start = [9 8];
goal = [3 1];
nexts = [
10 8;
7 9];
%}

figure;

test_kind = 3; % 1 = right tailed, 2 = left tailed, 3 = two-tailed


ms = [];
sems = [];
for t = 1:length(start)
    %which = s == start(t) & g == goal(t) & ord == ordinal(t); % ord is usually just 1
    which = s == start(t) & g == goal(t) & t_id == index(t); % for exp_v2_3 -- it compares trial indices proper
    move = dir(which);
    m = nexts(t,:);
    c1 = sum(move == m(1)); % count 1
    c2 = sum(move == m(2)); % count 2
    n = sum(which);

    switch test_kind 
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
    plot([0 3], [y(1) y(1)], '--', 'Color', [0.5 0.5 0.5]);
    plot([0 3], [y(2) y(2)], '--', 'Color', [0.5 0.5 0.5]);
    hold off;
    xticklabels({num2str(m(1)), num2str(m(2))});
    title(sprintf('%d -> %d: p = %.3f (c1 = %d, n = %d)', start(t), goal(t), p, c1, n));

    fprintf('trial #%d (%d -> %d): p = %.3f (c1 = %d, n = %d)\n', index(t), start(t), goal(t), p, c1, n);
    %ylim([4 5]);

    %{
    if t == 1
        ylabel('state chunking')
    elseif t == 3
        ylabel('action chunking / S-A')
    end
    %}
end

figure;
hold on;
bar(ms);
errorbar(ms, sems, 'linestyle', 'none', 'color', 'black');
plot([0 6], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 0.7], '-', 'color', [0.5 0.5 0.5])
hold off;
ylabel('p(go to 5)');
xticks(1:6);
xticklabels(index);
xlabel('trial # (probe)');
title(sprintf('human N = %d', length(data)));



%{
% do chunkers have shorter paths?
%
chunkers = s_id(ord == 1 & s == 6 & dir == 5);
nonchunkers = s_id(ord == 1 & s == 6 & dir == 7);
cl = subj_len(chunkers);
nl = subj_len(nonchunkers);
[h, p, ci, stats] = ttest2(cl, nl);
%}


clear all;


figure('pos', [10 500 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

modelfile = 'model_exp_v2_3_circ_alpha=1.0000_nsamples=10000_div_eps=0.6000_last.mat';
datafile = 'analyze_exp_v2_3.mat';


% A: graph
%
subplot(3,5,1);

load(modelfile);

H = chosen_H{1,end};
D = D_full(1);

c = [1 3 3 2 2 4 5 5 5 5];
[h, xs, ys] = plot_unlearn_graph(H, D, c);
%labelnode(h, 1:D.G.N, 1:D.G.N);
for i = 1:D.G.N
    text(h.XData(i) , h.YData(i) + 0.01, num2str(i), 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
set(gca, 'xlim', [-1 3]);
xlabel('stage 1');



subplot(3,5,2);

c = [3 3 3 2 2 2 5 5 5 5];
[h, xs, ys] = plot_unlearn_graph(H, D, c);
%labelnode(h, 1:D.G.N, 1:D.G.N);
for i = 1:D.G.N
    text(h.XData(i) , h.YData(i) + 0.01, num2str(i), 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
set(gca, 'xlim', [-1 3]);
xlabel('stage 2');

title('Experimental Design', 'fontsize', fontsize);


h = subplot(6,5,3);
pos = get(h, 'position');
pos(1) = pos(1) * 1.02;
pos(2) = pos(2) * 1.02;
pos(3) = pos(3) * 1.2;
pos(4) = pos(4) * 1.2;
subplot(6,5,3, 'position', pos);

PICpng = imread('exp/images/new_trial_crop.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  



h = subplot(6,5,8);
pos = get(h, 'position');
pos(1) = pos(1) * 1.02;
pos(2) = pos(2) * 1.02;
pos(3) = pos(3) * 1.2;
pos(4) = pos(4) * 1.2;
subplot(6,5,8, 'position', pos);
PICpng = imread('unlearn_crop.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  




h = subplot(3,3,3);
pos = get(h, 'position');
pos(1) = pos(1) * 0.93;
pos(2) = pos(2) * 0.98;
pos(3) = pos(3) * 1.3;
pos(4) = pos(4) * 1.3;
subplot(3,3, 3, 'position', pos);
PICpng = imread('unlearn_trials.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  




% B: Data
%


subplot(3,2,3);

load(datafile);

% swap action to be consistent with other plots
ms = 1 - ms;

stages = {[1:3], [4:6]}; % split in two

hold on;
for i = 1:2
    bar(stages{i}, ms(stages{i}));
end
for i = 1:2
    errorbar(stages{i}, ms(stages{i}), sems(stages{i}), 'linestyle', 'none', 'color', 'black');
end
plot([0 7], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 1.0], '-', 'color', [0.5 0.5 0.5])
xlim([0 7]);
ylim([0 0.9]);
hold off;
assert(nexts(1,1) == 7);
ylabel('P(action 6 \rightarrow 5)');
xticks(1:6);
xticklabels(1:6);
xlabel('probe trial');
legend({'stage 1', 'stage 2'}, 'location', 'northwest');

hold off;


title('Data', 'fontsize', fontsize);


% C: Model

load(modelfile);

subplot(3,2,4);

% swap action to be consistent with other plots
ms = 1 - ms;

hold on;
for i = 1:2
    bar(stages{i}, ms(stages{i}));
end
for i = 1:2
    errorbar(stages{i}, ms(stages{i}), sems(stages{i}), 'linestyle', 'none', 'color', 'black');
end
plot([0 7], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 1.0], '-', 'color', [0.5 0.5 0.5])
xlim([0 7]);
ylim([0 0.9]);
hold off;
assert(nexts(1,1) == 7);
xticks(1:6);
xticklabels(1:6);
xlabel('probe trial');
legend({'stage 1', 'stage 2'}, 'location', 'northwest');



title('Model', 'fontsize', fontsize);


% D: Hierarchies

%{
idx = [1 2 3 4 5 6 1 2 3 4 5 6];
subj = [1 1 1 1 1 1 2 2 2 2 2 2];
for i = 1:12
    subplot(6,6, 6*4 + i);

    s = subj(i);

    H = chosen_H{s, idx(i)};
    D = D_full(s);
    h = plot_H(H, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-2 4]);
    h.MarkerSize = 6;

    if i == 1
        ylabel('subject 1');
    end
    if i == 7
        ylabel('subject 2');
        xlabel('probe #1');
    end
    if i >= 8
        xlabel(sprintf('probe #%d', i - 6));
    end
end
%}



ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.09, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.09, 0.67, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.53, 0.67, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.09, 0.35, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');



% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
%set(h, 'PaperOrientation', 'landscape');
%print('figures/unlearn.pdf', '-dpdf');



%
%             stats
%


fprintf('\n\n --------------- DATA -----------------\n\n');

load(datafile);

% 3rd probe trial
%
which = t_id == index(3);
move = dir(which);
c1 = sum(move == nexts(3,1)); % count 1
c2 = sum(move ~= nexts(3,1)); % count 2
n = sum(which);
p = 2 * binocdf(min(c1,c2), n, 0.5);

fprintf('probe trial #3: %d out of %d participants, $p = %.2f$, two-tailed binomial test\n', c1, n, p);


% compare 3rd probe and 4th probe trial
%
which_g1 = t_id == index(3);
which_g2 = t_id == index(4);
assert(nexts(3,2) == nexts(4,2));  % to be consistent w/ other experiments

[h, p, ci, stats] = ttest2(dir(which_g1) == nexts(3,2), dir(which_g2) == nexts(4,2));

fprintf('probe trial #3 vs. #4: is there a difference? two-sample t-test: t(%d) = %.4f, p = %.4f\n', stats.df, stats.tstat, p);


% stats for LEARNING -- is there a ramp?
%
assert(nexts(3,2) == 5);
assert(nexts(4,2) == 5);
direction = dir == 5; % to be consistent w/ other experiments
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
fprintf('Learning: is there a ramp on probe trials 1..3? slope = %.2f, F(%d,%d) = %.2f, p = %.2f, mixed effects logistic regression with probe trials 1-3\n', H * beta, DF1, DF2, F, p);



% stats for UNLEARNing: trials 4..6
%
assert(nexts(3,2) == 5);
assert(nexts(4,2) == 5);
direction = dir == 5; % to be consistent w/ other experiments
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
assert(nexts(3,2) == 5);
assert(nexts(4,2) == 5);
direction = dir == 5; % to be consistent w/ other experiments

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
fprintf('Unlearning: is probe trial 4 below probe trial 3? slope = %.2f, F(%d,%d) = %.2f, p = %.10f, mixed effects logistic regression with probe trials 3-4\n', H * beta, DF1, DF2, F, p);




fprintf('\n\n --------------- MODEL -----------------\n\n');

load(modelfile);

% 3rd probe trial
%
move = mv(:,3);
n = length(move);
c1 = sum(move); %  count 1
c2 = sum(1 - move); % count 2
p = 2 * binocdf(min(c1,c2), n, 0.5);

fprintf('probe trial #3: %d out of %d simulated participants, $p = %e$, two-tailed binomial test\n', c1, n, p);

% compare 3rd probe and 4th probe trial
%
[h, p, ci, stats] = ttest2(mv(:,3), mv(:,4));

fprintf('probe trial #3 vs. #4: is there a difference? two-sample t-test: t(%d) = %.4f, p = %e\n', stats.df, stats.tstat, p);


% stats for LEARNING -- is there a ramp?
%
direction = 1 - mv(:,1:3); % to be consistent
direction = direction(:);
s = 1:length(D);
subject = repmat(s', [1 3]);
subject = subject(:);
trial_idx = repmat(1:3, [length(D), 1]);
trial_idx = trial_idx(:);
% create table
tbl = table(direction, trial_idx, subject);

formula = 'direction ~ 1 + trial_idx + (1 + trial_idx | subject)';
result1 = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result1);

H = [0 1];
[p, F, DF1, DF2] = coefTest(result1, H);
fprintf('Learning: is there a ramp on probe trials 1..3? slope = %.2f, F(%d,%d) = %.2f, p = %e, mixed effects logistic regression with probe trials 1-3\n', H * beta, DF1, DF2, F, p);


% stats for UNLEARNing: trials 4..6
%
direction = 1 - mv(:,4:6); % to be consistent
direction = direction(:);
s = 1:length(D);
subject = repmat(s', [1 3]);
subject = subject(:);
trial_idx = repmat(4:6, [length(D), 1]);
trial_idx = trial_idx(:);
% create table
tbl = table(direction, trial_idx, subject);

formula = 'direction ~ 1 + trial_idx + (1 + trial_idx | subject)';
result2 = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result2);

H = [1 0];
[p, F, DF1, DF2] = coefTest(result2, H); 



fprintf('Unlearning: are probe trials 4..6 below 0.5? intercept = %f, p = %e, F(%d,%d) = %f\n', H * beta, p, DF1, DF2, F);


% stats for UNLEARNing: trials 3 and 4
%
direction = 1 - mv(:,3:4); % to be consistent
direction = direction(:);
s = 1:length(D);
subject = repmat(s', [1 2]);
subject = subject(:);
trial_idx = repmat(3:4, [length(D), 1]);
trial_idx = trial_idx(:);
% create table
tbl = table(direction, trial_idx, subject);

formula = 'direction ~ 1 + trial_idx + (1 + trial_idx | subject)';
result3 = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result3);

H = [0 1];
[p, F, DF1, DF2] = coefTest(result3, H);
fprintf('Unlearning: is probe trial 4 below probe trial 3? slope = %.2f, F(%d,%d) = %.2f, p = %e, mixed effects logistic regression with probe trials 3-4\n', H * beta, DF1, DF2, F, p);

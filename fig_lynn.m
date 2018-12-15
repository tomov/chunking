clear all;
load('lynn_alpha=2_nsamples=1000.mat');

figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

% A: graph
%
subplot(2,3,1);

plot_lynn_graph(H, D);
title('Experimental Design', 'fontsize', fontsize);

hold on;
h = zeros(4,1);

h(1) = plot(NaN, NaN, 'o', 'color', [0.1 0.1 0.1], 'markerfacecolor', [255 255 255]/255);
h(2) = plot(NaN, NaN, 'o', 'color', [0.1 0.1 0.1], 'markerfacecolor', [29 117 0]/255);
h(3) = plot(NaN, NaN, 'o', 'color', [0.1 0.1 0.1], 'markerfacecolor', [0 118 186]/255);
h(4) = plot(NaN, NaN, 'o', 'color', [0.1 0.1 0.1], 'markerfacecolor', [181 23 0]/255);

hold off;
legend(h, {'0 (current node)', '1 (no violation)', '2 (short violation)', '3,4 (long violation)'}, 'Position', [0.06 0.82 0.05 0.1]);



% B: Data
%

subplot(2,3,2);

p_short = p_cross(:,2) - p_cross(:,1);
p_long = p_cross(:,3) - p_cross(:,1);

m1 = 38;
se1 = 7;
m2 = 63;
se2 = 6;

bar([m1 m2]);
hold on;
errorbar([m1 m2], [se1 se2], 'linestyle', 'none', 'color', 'black');
hold off;

%set(gca, 'xlim', [0 2]);
%set(gca, 'ylim', [0 1]);
%set(gca, 'ytick', [0 0.5 1]);
xticklabels({'short', 'long'});
h = xlabel('violations');
pos = get(h,'Position');
pos(2) = pos(2) * 0.7;
set(h, 'Position', pos);
ylabel('Change in RT (ms)');
hold off;

title('Data', 'fontsize', fontsize);


% C: Model

subplot(2,3,3);

p_short = p_cross(:,2) - p_cross(:,1);
p_long = p_cross(:,3) - p_cross(:,1);

m1 = mean(p_short);
se1 = sem(p_short);
m2 = mean(p_long);
se2 = sem(p_long);

bar([m1 m2]);
hold on;
errorbar([m1 m2], [se1 se2], 'linestyle', 'none', 'color', 'black');
hold off;

%set(gca, 'xlim', [0 2]);
%set(gca, 'ylim', [0 1]);
%set(gca, 'ytick', [0 0.5 1]);
%xticklabels({'P(fewer boundaries)'});
xticklabels({'short', 'long'});
h = xlabel('violations');
pos = get(h,'Position');
pos(2) = pos(2) * 0.7;
set(h, 'Position', pos);
ylabel('P(different cluster)');
hold off;

title('Model', 'fontsize', fontsize);


% D: Hierarchies

for s = 1:12
    subplot(4,6, 12 + s);
    h = plot_H(map_H{s}, D);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    h.MarkerSize = 6;

    if s == 3
        title(['                            ', 'Example hierarchies'], 'fontsize', fontsize);
    end
end


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.10, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.36, 0.96, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.64, 0.96, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.10, 0.52, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');


% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/lynn.pdf', '-dpdf');


% stats
%
[h, p, ci, stats] = ttest(p_short);
fprintf('short violations: t(%d) = %.2f, p = %e (one sample two-tailed t-test against 0)\n', stats.df, stats.tstat, p);

[h, p, ci, stats] = ttest(p_long);
fprintf('long violations: t(%d) = %.2f, p = %e (one sample two-tailed t-test against 0)\n', stats.df, stats.tstat, p);

[h, p, ci, stats] = ttest2(p_long, p_short);
fprintf('long vs. short violations: t(%d) = %.2f, p = %e (two sample two-tailed t-test)\n', stats.df, stats.tstat, p);

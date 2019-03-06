function fig_solway2(filename, do_save)

if ~exist('filename', 'var') || isempty(filename)
    %load('schapiro_alpha=2.mat'); % nsamples = 1000
    %load('schapiro_alpha=2_nsamples=20.mat');

    %load('schapiro_alpha=2_nsamples=100.mat'); % <-- preprint
    load schapiro_N=30_alpha=2.0000_nsamples=1000.mat  % <-- sample_c
else
    load(filename);
end

if ~exist('do_save', 'var')
    do_save = false;
end

figure('pos', [10 1200 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

% A: graph
%
subplot(2,3,1);

r1 = 6;
r2 = 3;
cx = cos([0:2]*2*pi/3 + pi/2) * r1;
cy = sin([0:2]*2*pi/3 + pi/2) * r1;
x(1:5) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(1);
y(1:5) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(1);
x(6:10) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(2);
y(6:10) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(2);
x(11:15) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(3);
y(11:15) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(3);


H.c = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3];

h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    highlight(h, i, 'nodecolor', [0.1 0.1 0.1], 'markersize', 12);
end
hold on;

h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    if H.c(i) == 1
        highlight(h, i, 'nodecolor', [113 84 159]/255);
    elseif H.c(i) == 2
        highlight(h, i, 'nodecolor', [225 113 39]/255);
    else
        highlight(h, i, 'nodecolor', [5 145 72]/255);
    end

    for j = i+1:D.G.N
        if D.G.E(i,j)
            highlight(h, i, j, 'edgecolor', [0.1 0.1 0.1], 'linewidth', 2);
        end
    end
end

set(gca, 'xtick', []);
set(gca, 'ytick', []);
title('Experimental Design', 'fontsize', fontsize);


% B: Data

subplot(2,3,2);

m = [0.305 0.245;
     0.3   0.246];
se = [0.025; 0.026;
      0.025; 0.024];

h = bar(m);
xs = [h(1).XData + h(1).XOffset, ...
      h(2).XData + h(2).XOffset];
h(1).FaceColor = 'flat';
h(2).FaceColor = 'flat';
h(1).CData(1,:) = [0.3 0.3 0.3];
h(1).CData(2,:) = h(1).CData(1,:);
h(2).CData(1,:) = [0.6 0.6 0.6];
h(2).CData(2,:) = h(2).CData(1,:);
hold on;
m = m(:); se = se(:); 
errorbar(xs, m, se, 'linestyle', 'none', 'color', 'black');
hold off;

xticklabels({'All trials', 'Hamiltonian paths'});
set(gca,'FontSize', axisfontsize);
ylabel('Probability of parse', 'fontsize', axisfontsize);
title('Data', 'fontsize', fontsize);


% C: Model

subplot(2,3,3);

m = [mean(comm_p) mean(other_p);
     mean(comm_p_hamil) mean(other_p_hamil)];
se = [sem(comm_p) sem(other_p);
      sem(comm_p_hamil) sem(other_p_hamil)];

h = bar(m);
xs = [h(1).XData + h(1).XOffset, ...
      h(2).XData + h(2).XOffset];
h(1).FaceColor = 'flat';
h(2).FaceColor = 'flat';
h(1).CData(1,:) = [0.3 0.3 0.3];
h(1).CData(2,:) = h(1).CData(1,:);
h(2).CData(1,:) = [0.6 0.6 0.6];
h(2).CData(2,:) = h(2).CData(1,:);
hold on;
m = m(:); se = se(:); 
errorbar(xs, m, se, 'linestyle', 'none', 'color', 'black');
hold off;

xticklabels({'All trials', 'Hamiltonian paths'});
lgd = legend({['Community' char(10) 'transition parse'], 'Other parse'});
lgd.Position(1) = 0.86;
title('Model', 'fontsize', fontsize);
set(gca,'FontSize', axisfontsize);

% D: Hierarchies

for s = 1:12
    subplot(4,6, 12 + s);
    h = plot_H(map_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);

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
if do_save
    h = gcf;
    %set(h, 'PaperPositionMode', 'auto');
    set(h, 'PaperOrientation', 'landscape');
    print('figures/schapiro.pdf', '-dpdf');
end


% stats
%

[~, p, ci, stats] = ttest2(comm_p, other_p);
fprintf('random walks: t(%d) = %.2f, p = %e (two sample two-tailed t-test)\n', stats.df, stats.tstat, p);

[~, p, ci, stats] = ttest2(comm_p_hamil, other_p_hamil);
fprintf('hamiltonians: t(%d) = %.2f, p = %e (two sample two-tailed t-test)\n', stats.df, stats.tstat, p);



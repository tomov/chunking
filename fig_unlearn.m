clear all;


figure('pos', [10 500 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

modelfile = 'model_exp_v2_3_circ_alpha=1.0000_nsamples=10000_div_eps=0.6000_last.mat';

% A: graph
%
subplot(3,6,1);

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
xlabel('bad');



subplot(3,6,2);

[h, xs, ys] = plot_unlearn_graph(H, D);
%labelnode(h, 1:D.G.N, 1:D.G.N);
for i = 1:D.G.N
    text(h.XData(i) , h.YData(i) + 0.01, num2str(i), 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
set(gca, 'xlim', [-1 3]);
xlabel('good');

title('Experimental Design', 'fontsize', fontsize);


h = subplot(6,6,3);
pos = get(h, 'position');
pos(1) = pos(1) * 0.96;
pos(2) = pos(2) * 0.98;
pos(3) = pos(3) * 1.2;
pos(4) = pos(4) * 1.2;
subplot(6,6,3, 'position', pos);

PICpng = imread('exp/images/new_trial_crop.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  



h = subplot(6,6,9);
pos = get(h, 'position');
pos(1) = pos(1) * 0.95;
pos(2) = pos(2) * 1.0;
pos(3) = pos(3) * 1.2;
pos(4) = pos(4) * 1.2;
subplot(6,6,9, 'position', pos);
PICpng = imread('unlearn_crop.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  




h = subplot(3,2,2);
pos = get(h, 'position');
pos(1) = pos(1) * 0.90;
pos(2) = pos(2) * 0.98;
pos(3) = pos(3) * 1.3;
pos(4) = pos(4) * 1.3;
subplot(3,2, 2, 'position', pos);
PICpng = imread('subway9_map_trials.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  




% B: Data
%


subplot(3,2,3);

load('analyze_exp_v2_3.mat');

% swap action to be consistent with other plots
ms = 1 - ms;

hold on;
bar(ms);
errorbar(ms, sems, 'linestyle', 'none', 'color', 'black');
plot([0 6], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 0.8], '-', 'color', [0.5 0.5 0.5])
hold off;
assert(nexts(1,1) == 7);
ylabel('P(action 6 \rightarrow 5)');
xticks(1:6);
xticklabels(1:6);
xlabel('probe trial');

hold off;


title('Data', 'fontsize', fontsize);


% C: Model

load(modelfile);

subplot(3,2,4);

% swap action to be consistent with other plots
ms = 1 - ms;

hold on;
bar(ms);
errorbar(ms, sems, 'linestyle', 'none', 'color', 'black');
plot([0 6], [0.5 0.5], '--', 'color', [0.5 0.5 0.5])
plot([3.5 3.5], [0 0.8], '-', 'color', [0.5 0.5 0.5])
hold off;
assert(nexts(1,1) == 7);
xticks(1:6);
xticklabels(1:6);
xlabel('probe trial');



title('Model', 'fontsize', fontsize);


% D: Hierarchies

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



ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.09, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.09, 0.67, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.53, 0.67, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.09, 0.35, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');



% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/unlearn.pdf', '-dpdf');




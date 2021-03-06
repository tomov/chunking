clear all;


figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

%modelfile = 'trash/model_all_data_100samples_MAP_5alpha.mat';
modelfile = 'model_all_data_40samples_MAP_2alpha.mat';

ii = 4;
jj = 1:3;

% A: graph
%
subplot(2,5,1);

load(modelfile);

H = pl(ii).H{jj(1)}(1,1);
D = pl(ii).D{jj(1)}(1,1);
[h, xs, ys] = plot_subway9_graph(H, D);
labelnode(h, 1:D.G.N, 1:D.G.N);





% B: Data
%


subplot(2,5,4);

load('analyze_all_data.mat');

i = ii;
m = pl(i).m ./ pl(i).n;
ci = pl(i).ci ./ pl(i).n;
for j = 1:length(pl(i).m)
    se(j) = std([ones(1,pl(i).m(j)) zeros(1,pl(i).n(j) - pl(i).m(j))]) / sqrt(pl(i).n(j));
end

hold on;
col = [];
for j = jj
    h = bar(j, m(j));
    if pl(i).fake(j)
        set(h, 'facealpha', 0.5);
    end
    if isempty(col)
        col = get(h, 'facecolor');
    else
        set(h, 'facecolor', col);
    end

    %text(j, 1, sprintf('p = %.3f\nN = %d', pl(i).p(j), pl(i).n(j)));
    %text(j - 0.3, 0.9 - j * 0.1, sprintf('p = %.3f', pl(i).p(j)));
end
errorbar(m(jj), se(jj), 'linestyle', 'none', 'color', 'black');
%errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
line([0 length(jj)+1], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 length(jj)+1]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
xticklabels({'P(fewer boundaries)'});
ylabel('fraction of participants');

hold off;



title('Data', 'fontsize', fontsize);


% C: Model

load(modelfile);

subplot(2,5,5);

i = ii;
m = pl(i).m ./ pl(i).n;
ci = pl(i).ci ./ pl(i).n;
for j = 1:length(pl(i).m)
    se(j) = std([ones(1,pl(i).m(j)) zeros(1,pl(i).n(j) - pl(i).m(j))]) / sqrt(pl(i).n(j));
end

hold on;
col = [];
for j = jj
    h = bar(j, m(j));
    if pl(i).fake(j)
        set(h, 'facealpha', 0.5);
    end
    if isempty(col)
        col = get(h, 'facecolor');
    else
        set(h, 'facecolor', col);
    end

    %text(j, 1, sprintf('p = %.3f\nN = %d', pl(i).p(j), pl(i).n(j)));
   % text(j - 0.3, 0.9 - j * 0.1, sprintf('p = %.3f', pl(i).p(j)));
end
errorbar(m(jj), se(jj), 'linestyle', 'none', 'color', 'black');
%errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
line([0 length(jj)+1], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 length(jj)+1]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
xticklabels({'P(fewer boundaries)'});
ylabel('fraction of participants');



title('Model', 'fontsize', fontsize);


% D: Hierarchies

j = [1 1 1 1 1 1 2 2 2 2 2 2];
for s = 1:12
    subplot(4,6, 12 + s);

    H = pl(i).H{j(s)}(s,:);
    D = pl(i).D{j(s)}(s);
    P = pl(i).P{j(s)}(s,:);
    [~,I] = max(P); % MAP H
    H = H(I);
    map_H{s} = H;
    h = plot_H(map_H{s}, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-2 4]);
    h.MarkerSize = 6;

    if s == 3
        title('Example hierarchies', 'fontsize', fontsize);
    end
    if s == 1
        ylabel('bad');
    end
    if s == 7
        ylabel('control');
    end
end

% A: design; it's here b/c otherwise it disappears...

h = subplot(2,5,2);
pos = get(h, 'position');
pos(1) = pos(1) * 0.9;
pos(2) = pos(2) * 0.6;
pos(3) = pos(3) * 2.4;
pos(4) = pos(4) * 2.4;
subplot(2,5, 2, 'position', pos);
PICpng = imread('subway9_trials.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  

title('Experimental Design', 'fontsize', fontsize);


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.10, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.59, 0.96, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.75, 0.96, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.10, 0.52, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');



% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/subway9.pdf', '-dpdf');




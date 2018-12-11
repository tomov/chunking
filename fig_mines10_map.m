clear all;


figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

ii = 5;
jj = 1;

% A: graph
%
subplot(2,6,1);

load('model_all_data_10samples_MAP_5alpha.mat');

H = pl(1).H{jj}(1,1);
D = pl(1).D{jj}(1,1);
[h, xs, ys] = plot_subway10_graph(H, D);
labelnode(h, 1:D.G.N, 1:D.G.N);




h = subplot(2,6,2);
pos = get(h, 'position');
pos(1) = pos(1) * 1.0;
pos(2) = pos(2) * 1.0;
pos(3) = pos(3) * 1.0;
pos(4) = pos(4) * 1.0;
subplot(2,6, 2, 'position', pos);

PICpng = imread('mines10_free_crop.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  

xlabel('free choice');





h = subplot(2,6,3);
pos = get(h, 'position');
pos(1) = pos(1) * 1.0;
pos(2) = pos(2) * 1.0;
pos(3) = pos(3) * 1.0;
pos(4) = pos(4) * 1.0;
subplot(2,6, 3, 'position', pos);

PICpng = imread('mines10_forced_crop.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  
xlabel('forced choice');

title('Experimental Design', 'fontsize', fontsize);



h = subplot(2,6,4);
pos = get(h, 'position');
pos(1) = pos(1) * 1.0;
pos(2) = pos(2) * 1.0;
pos(3) = pos(3) * 1.0;
pos(4) = pos(4) * 1.0;
subplot(2,6, 4, 'position', pos);

PICpng = imread('mines10_trials.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  


% B: Data
%


subplot(2,6,5);

load('analyze_all_data.mat');

i = ii;
m = pl(i).m ./ pl(i).n;
ci = pl(i).ci ./ pl(i).n;
for j = 1:length(pl(i).m)
    se(j) = std([ones(1,pl(i).m(j)) zeros(1,pl(i).n(j) - pl(i).m(j))]) / sqrt(pl(i).n(j));
end

hold on;
j = jj;
h = bar(j, m(j));
errorbar(m(j), se(j), 'linestyle', 'none', 'color', 'black');
%errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
xticklabels({'P(fewer boundaries)'});
ylabel('fraction of participants');

hold off;



title('Data', 'fontsize', fontsize);


% C: Model

%load('model_all_data_10samples_MAP_5alpha.mat');

subplot(2,6,6);

i = ii;

m = 0.8; % TODO what is the fraction of simulations that go to state 5?
se = 0.04; % TODO what is the standard error of the mean for it?
%m = pl(i).m ./ pl(i).n;
%ci = pl(i).ci ./ pl(i).n;
%for j = 1:length(pl(i).m)
%    se(j) = std([ones(1,pl(i).m(j)) zeros(1,pl(i).n(j) - pl(i).m(j))]) / sqrt(pl(i).n(j));
%end

hold on;
j = jj;
h = bar(j, m(j));
errorbar(m(j), se(j), 'linestyle', 'none', 'color', 'black');
%errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
xticklabels({'P(fewer boundaries)'});
hold off;


title('Model', 'fontsize', fontsize);


% D: Hierarchies

%{
TODO
for s = 1:12
    subplot(4,6, 12 + s);

    H = pl(i).H{j}(s,:);
    D = pl(i).D{j}(s);
    P = pl(i).P{j}(s,:);
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
end
%}


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.10, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.64, 0.96, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.80, 0.96, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.10, 0.52, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');


% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/mines10_map.pdf', '-dpdf');




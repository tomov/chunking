% figure for experiment 6
%
clear all;


figure('pos', [2000 500 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 14;

%modelfile = 'mines10_alpha=2_nsamples=1000.mat'; % <--- preprint
%modelfile = 'mines10_alpha=1_nsamples=10000_last.mat'; <-- sample_c
%modelfile = 'mines10_alpha=1.0000_nsamples=10000_last.mat'; % sample_c
modelfile = 'mines10_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat';

ii = 5;
jj = 1;

% A: graph
%
subplot(2,6,1);

load(modelfile);
[h, xs, ys] = plot_subway10_graph(H, D);
%labelnode(h, 1:D.G.N, 1:D.G.N);
for i = 1:D.G.N
    text(h.XData(i) , h.YData(i) + 0.01, num2str(i), 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end


title(['                                                                                              ','Experimental Design'], 'fontsize', fontsize);


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




h = subplot(2,6,4);
pos = get(h, 'position');
pos(1) = pos(1) * 0.96;
pos(2) = pos(2) * 0.97;
pos(3) = pos(3) * 1.2;
pos(4) = pos(4) * 1.2;
subplot(2,6, 4, 'position', pos);

PICpng = imread('mines10_trials.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  


% B: Data
%


subplot(2,6,5);

load('analyze_Exp_1_thru_5.mat');

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
%xticklabels({'P(fewer boundaries)'});
%ylabel('fraction of participants');
xticklabels({'fraction participants'});
ylabel('P(action 6 \rightarrow 5)');

hold off;

assert(pl(i).tests(j) == 3);
%fprintf('two-tailed binomial test c1 = %d, n = %d, p = %e\n', pl(i).m(j), pl(i).n(j), pl(i).p(j));
fprintf('%d out of %d participants, $p = %.2f$, two-tailed binomial test\n', pl(i).m(j), pl(i).n(j), pl(i).p(j));


title('Data', 'fontsize', fontsize);


% C: Model

load(modelfile);

subplot(2,6,6);


hold on;
h = bar(m);
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
%errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%xticklabels({'P(fewer boundaries)'});
xticklabels({'fraction simulations'});
hold off;

%fprintf('two-tailed binomial test c1 = %d, n = %d, p = %e\n', c1, n, p);
fprintf('%d out of %d simulated participants, $p = %.2f$, two-tailed binomial test\n', pl(i).m(j), pl(i).n(j), pl(i).p(j));

title('Model', 'fontsize', fontsize);


% D: Hierarchies

%{
for s = 1:12
    subplot(4,6, 12 + s);

    h = plot_H(chosen_H{s}, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-2 4]);
    h.MarkerSize = 6;

    if s == 3
        title(['                            ', 'Example hierarchies'], 'fontsize', fontsize);
    end
end
%}


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.10, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.64, 0.96, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.80, 0.96, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.10, 0.52, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');


% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/mines10_map.pdf', '-dpdf');




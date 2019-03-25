% figure for experiment 7
%
clear all;


%HUMAN DATA
%
%<First Plot: graphs 1, 4, 2, 5>
%means = [34/40, 34/40, 19/40, 10/40]
%errors = [0.05717718748968655, 0.05717718748968655, 0.07996393417804533, 0.06933752452815363]
%
%<Second Plot: graphs 6, 3, 7>
%means = [37/40, 31/40, 6/40]
%errors = [0.04217636961434867, 0.06686668711812965, 0.05717718748968655]
%
%
%MODEL
%
%<First Plot: graphs 1, 4, 2, 5>
%means = [20/40, 30/40, 3/40, 5/40]
%errors = [0.08006407690254357, 0.06933752452815363, 0.04217636961434867, 0.05295740910852021]
%
%<Second Plot: graphs 6, 3, 7>
%means = [33/40, 27/40, 11/40]
%errors = [0.06084343084444759, 0.07499999999999998, 0.07149950690165274]



figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;


% A: graph
%


h = subplot(2,1,1);
pos = get(h, 'position');
pos(1) = pos(1) * 0.5;
pos(2) = pos(2) * 1.0;
pos(3) = pos(3) * 1.2;
pos(4) = pos(4) * 1.2;
subplot(2,1, 1, 'position', pos);


PICpng = imread('active_reordered.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  

title('Experimental Design');


% B: Data
%


human_data = [34, 19, 31, 34, 10, 37, 6];
graphs_with_two_edges = [1 4 2 5];
graphs_with_three_edges = [6 3 7];

% we show the graphs in a weird order
c_two = human_data(graphs_with_two_edges);
c_three = human_data(graphs_with_three_edges);
n = 40;

m_two = c_two / n;
m_three = c_three / n;

for i = 1:length(c_two)
    c = c_two(i);
    se_two(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_two(i) = binocdf(min(c,n-c),n,0.5) + 1 - binocdf(max(c,n-c),n,0.5);
end
for i = 1:length(c_three)
    c = c_three(i);
    se_three(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_three(i) = 1 - binocdf(c,n,1/3);
end

m_data = [m_two m_three];
p_data = [p_two p_three];
p_data_corr = 1 - (1 - p_data) .^ length(p_data);

disp('p-values for data (uncorrected)');
p_data
disp('p-values for data (Bonferroni corrected)');
p_data_corr



subplot(2,4,5);


%hold on;
%m = [m_two m_three];
%se = [se_two se_three];
%x = [1:4 6:8];
%bar(x, m);
%errorbar(x, m, se, 'linestyle', 'none', 'color', 'black');
%line([0 5], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%line([5 9], [1/3 1/3], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%set(gca, 'xlim', [0 9]);
%set(gca, 'ylim', [0 1]);
%%set(gca, 'ytick', [0 0.5 1]);
%set(gca, 'xtick', [1 2 3 4 6 7 8]);
%xticklabels([1 2 3 4 5 6 7]);
%%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%%xticklabels({'P(fewer boundaries)'});
%ylabel('fraction of participants');
%xlabel('graph');



hold on;
bar(m_two);
errorbar(m_two, se_two, 'linestyle', 'none', 'color', 'black');
line([0 5], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 5]);
set(gca, 'ylim', [0 1]);
%set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1 2 3 4]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%xticklabels({'P(fewer boundaries)'});
ylabel('P(edge A)');
xlabel('graph');

hold off;

title(['                                          ', 'Data'], 'fontsize', fontsize);

subplot(2,4,6);

hold on;
bar(m_three);
errorbar(m_three, se_three, 'linestyle', 'none', 'color', 'black');
line([0 4], [1/3 1/3], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 4]);
set(gca, 'ylim', [0 1]);
%set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1 2 3]);
xticklabels([5 6 7]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%xticklabels({'P(fewer boundaries)'});
ylabel('P(edge B)');
xlabel('graph');

hold off;








% C: Model

% from Wanqian:
% model is in https://github.com/Wanqianxn/hierarchy-edge-unveiling
%
% Human     {34, 19, 31, 34, 10, 37, 06} # human data from last year
% 1.0e1.0   {36, 00, 29, 37, 00, 31, 12} # model data without tiebreaks, epsilon=1.0 (no exploring)
% 1.0e0.8   {33, 03, 26, 34, 05, 30, 11} # model data without tiebreaks, epsilon=0.8
% 1.0e0.6   {32, 09, 29, 33, 07, 24, 10} # model data without tiebreaks, epsilon=0.6
% 1.0te1.0  {20, 00, 29, 26, 00, 31, 12} # model data with tiebreaks, epsilon=1.0 (no exploring)
% 1.0te0.8  {15, 05, 24, 28, 07, 27, 12} # model data with tiebreaks, epsilon=0.8
% 1.0te0.6  {17, 12, 22, 28, 09, 24, 16} # model data with tiebreaks, epsilon=0.6 % <-- this is it

h = init_hyperparams;
assert(abs(h.alpha - 1) < 1e-10);
assert(abs(h.eps - 0.6) < 1e-10);

model_data = [17, 12, 22, 28, 09, 24, 16]; % from Wan's simulations with alpha = 1, eps = 0.6


c_two = model_data(graphs_with_two_edges);
c_three = model_data(graphs_with_three_edges);
n = 40;

m_two = c_two / n;
m_three = c_three / n;

for i = 1:length(c_two)
    c = c_two(i);
    se_two(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_two(i) = binocdf(min(c,n-c),n,0.5) + 1 - binocdf(max(c,n-c),n,0.5);
end
for i = 1:length(c_three)
    c = c_three(i);
    se_three(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_three(i) = 1 - binocdf(c,n,1/3);
end


m_model = [m_two m_three];
p_model = [p_two p_three];
p_model_corr = 1 - (1 - p_model) .^ length(p_model);

disp('p-values for model (uncorrected)');
p_model
disp('p-values for model (Bonferroni corrected)');
p_model_corr


subplot(2,4,7);


%hold on;
%m = [m_two m_three];
%se = [se_two se_three];
%x = [1:4 6:8];
%bar(x, m);
%errorbar(x, m, se, 'linestyle', 'none', 'color', 'black');
%line([0 5], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%line([5 9], [1/3 1/3], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%set(gca, 'xlim', [0 9]);
%set(gca, 'ylim', [0 1]);
%%set(gca, 'ytick', [0 0.5 1]);
%set(gca, 'xtick', [1 2 3 4 6 7 8]);
%xticklabels([1 2 3 4 5 6 7]);
%%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%%xticklabels({'P(fewer boundaries)'});
%ylabel('fraction of participants');
%xlabel('graph');



hold on;
bar(m_two);
errorbar(m_two, se_two, 'linestyle', 'none', 'color', 'black');
line([0 5], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 5]);
set(gca, 'ylim', [0 1]);
%set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1 2 3 4]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%xticklabels({'P(fewer boundaries)'});
ylabel('P(edge A)');
xlabel('graph');

hold off;

title(['                                          ', 'Model'], 'fontsize', fontsize);


subplot(2,4,8);

hold on;
bar(m_three);
errorbar(m_three, se_three, 'linestyle', 'none', 'color', 'black');
line([0 4], [1/3 1/3], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [0.5 - ci(j) 0.5 - ci(j) 0.5 + ci(j) 0.5 + ci(j)], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 4]);
set(gca, 'ylim', [0 1]);
%set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1 2 3]);
xticklabels([5 6 7]);
%text(0.7, 0.9, sprintf('p = %.3f', pl(i).p(j)));
%xticklabels({'P(fewer boundaries)'});
ylabel('P(edge B)');
xlabel('graph');

hold off;






[r,p] = corr(m_data', m_model');
fprintf('Pearson r = %.4f, p = %.4f\n', r, p);





ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.10, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.10, 0.52, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.50, 0.52, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');


% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/active.pdf', '-dpdf');




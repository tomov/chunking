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


% B: Data
%


c_two = [34, 34, 19, 10];
c_three = [37, 31, 6];
n = 40;

m_two = c_two / n;
m_three = c_three / n;

for i = 1:length(c_two)
    c = c_two(i);
    se_two(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_two(i) = 1 - binocdf(min(c,n-c),n,0.5);
end
for i = 1:length(c_three)
    c = c_three(i);
    se_three(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_three(i) = 1 - binocdf(min(c,n-c),n,0.5);
end



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


title('Data', 'fontsize', fontsize);






% C: Model






c_two = [20 30 3 5];
c_three = [33 27 11];
n = 40;

m_two = c_two / n;
m_three = c_three / n;

for i = 1:length(c_two)
    c = c_two(i);
    se_two(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_two(i) = 1 - binocdf(min(c,n-c),n,0.5);
end
for i = 1:length(c_three)
    c = c_three(i);
    se_three(i) = std([ones(1,c) zeros(1,n - c)]) / sqrt(n);
    p_three(i) = 1 - binocdf(min(c,n-c),n,0.5);
end



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


title('Model', 'fontsize', fontsize);









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




clear all;


figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

% A: graph
%
subplot(2,3,1);

% TODO load from model output file 
D = init_D_from_txt('mines.txt');
h.alpha = 5;
H = init_H(D,h);
[h, xs, ys] = plot_mines_graph(H, D);
labelnode(h, 1:D.G.N, 1:D.G.N);





% B: Data
%


subplot(2,3,2);


% TODO from data file 
c = 17;
n = 20;
m = c / 20;
se = std([ones(1,c) zeros(1,n-c)]) / sqrt(n);
p = 1 - binocdf(c,n,0.5);

hold on;
h = bar(m);
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
text(0.7, 0.9, sprintf('p = %.3f', p));
xticklabels({'P(left)'});
ylabel('fraction of participants');

hold off;



title('Data', 'fontsize', fontsize);


% C: Model

subplot(2,3,3);


% TODO from model output 
c = 16;
n = 20;
m = c / 20;
se = std([ones(1,c) zeros(1,n-c)]) / sqrt(n);
p = 1 - binocdf(c,n,0.5);

fprintf('right-tailed binomial test c = %d, n = %d, p = %.4f\n', c, n, p);

hold on;
h = bar(m);
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
text(0.7, 0.9, sprintf('p = %.3f', p));
xticklabels({'P(left)'});

hold off;


title('Model', 'fontsize', fontsize);


% D: Hierarchies

% TODO load MAP hierarchies from model output
%{
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
text(0.36, 0.96, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.65, 0.96, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.10, 0.52, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');


% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('mines.pdf', '-dpdf');



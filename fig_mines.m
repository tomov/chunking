clear all;


figure('pos', [1000 500 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

%modelfile = 'mines_alpha=2_nsamples=1000.mat'; % <--- preprint
%modelfile = 'mines_alpha=1_nsamples=10000_last.mat'; <-- sample_c
%modelfile = 'mines_alpha=1.0000_nsamples=10000_last.mat'; % sample_c
modelfile = 'mines_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat'; % sample_c

% A: graph
%
subplot(2,3,1);

load(modelfile);
[h, xs, ys] = plot_mines_graph(H, D);
%labelnode(h, 1:D.G.N, 1:D.G.N);
for i = 1:D.G.N
    text(h.XData(i) , h.YData(i) + 0.01, num2str(i), 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end





% B: Data
%


subplot(2,3,2);



% from Agni & Samyu's paper
c1 = 24;
n = 32;
c2 = n - c1;
m = c1 / n;
se = std([ones(1,c1) zeros(1,c2)]) / sqrt(n);
p = 2 * binocdf(min(c1,c2),n,0.5);

fprintf('Data: two-tailed binomial test c = %d, n = %d, p = %.4f\n', c1, n, p);

hold on;
h = bar(m);
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
%text(0.7, 0.9, sprintf('p = %.3f', p));
%xticklabels({'P(state 3)'});
xticklabels({'fraction participants'});
%ylabel('fraction of participants');
ylabel('P(choose 3)');

hold off;


title('Data', 'fontsize', fontsize);



% C: Model

subplot(2,3,3);

load(modelfile);

fprintf('Model: two-tailed binomial test c1 = %d, n = %d, p = %.4f\n', c1, n, p);

hold on;
h = bar(m);
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
set(gca, 'xtick', [1]);
%text(0.7, 0.9, sprintf('p = %.3f', p));
%xticklabels({'P(state 3)'});
xticklabels({'fraction simulations'});

hold off;


title('Model', 'fontsize', fontsize);


% D: Hierarchies

for s = 1:12
    subplot(4,6, 12 + s);

    h = plot_H(chosen_H{s}, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'ylim', [-3 3]);
    h.MarkerSize = 6;

    if s == 3
        title(['                            ', 'Example hierarchies'], 'fontsize', fontsize);
    end
end


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
print('figures/mines.pdf', '-dpdf');



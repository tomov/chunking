% Solway experiment 1 visualization
clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

h.alpha = 5;

D = init_D_from_txt('solway1.txt');

for i = 1:10
    i
    [H, P] = sample(D, h, 1000);
    H_all{i} = H;
    P_all{i} = P;

    [~,I] = max(P);
    H_map(i) = H(I);
end

save('solway1_viz.mat');

load('solway1_viz.mat');

x = [-3 -1 -2 -3 -1 1 3 2 1 3];
y = [-3 -3 -2 -1 -1 1 1 2 3 3];
y = -y;

figure;

for i = 1:length(H_map)
    subplot(1, length(H_map), i);
    h = plot_H(H_map(i), D);
    set(h, 'XData', x);
    set(h, 'YData', y);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
end

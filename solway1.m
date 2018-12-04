% simulate experiment 1 from solway 2014

clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 40; % participants
h.alpha = 2;
nsamples = 100;

D = init_D_from_txt('solway1.txt');


for s = 1:N % for each simulated subject
    fprintf('subject %d\n', s);

    [H, P] = sample(D, h, nsamples);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);
    map_H{s} = H;

    H = populate_H(H, D); % fill up bridges

    b = zeros(1,D.G.N);
    for i = 1:length(H.cnt)
        for j = i+1:length(H.cnt)
            bridge = H.b{i,j};
            if ~isempty(bridge)
                b(bridge(1)) = 1;
                b(bridge(2)) = 1;
            end
        end
    end

    b = find(b);
    if isempty(b)
        % only 1 cluster -> all nodes are fine
        b = 1:D.G.N;
    end
    loc(s,:) = datasample(b, 3); % pick 3 at random (with replacement)
end


filename = sprintf('solway1_alpha=%d_nsamples=%d.mat', h.alpha, nsamples);
save(filename);

%load('solway1.mat');

x = [-3 -1 -2 -3 -1 1 3 2 1 3];
y = [-3 -3 -2 -1 -1 1 1 2 3 3];
y = -y;

figure;

H.c = [1 1 1 1 1 2 2 2 2 2];
h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
for i = 1:D.G.N
    highlight(h, i, 'NodeColor', [0.6 0.6 0.6], 'MarkerSize', 10 + sum(loc(:) == i));
end
set(gca, 'xtick', []);
set(gca, 'ytick', []);

hold on;
D.G.E(:) = 0;
h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
hold off;


c1 = sum(loc(:) == 5 | loc(:) == 6); % count 1
c2 = sum(loc(:) ~= 5 & loc(:) ~= 6); % count 2
n = c1 + c2;
p = 1 - binocdf(c1, n, 2/10);
y = binoinv([0.025 0.975], n, 2/10) / n;

figure;
m = c1/n;
se = std(loc(:) == 5 | loc(:) == 6) / sqrt(n);
bar(m);
hold on;
errorbar(m, se);
line([0 2], [2/10 2/10], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [y(1) y(1) y(2) y(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');
hold off;


fprintf('right-tailed binomial test n = %d, p = %.4f\n', n, p);

% simulate experiment 1 from solway 2014

function filename = solway1(N, h, nsamples, take_map)

rng default;

sem = @(x) std(x) / sqrt(length(x));

if ~exist('N', 'var') || isempty(N)
    N = 40; % participants
end
if ~exist('h', 'var')
    h = init_hyperparams;
end
if ~exist('nsamples', 'var')
    nsamples = 10000;
end
if ~exist('take_map', 'var')
    take_map = false;
end

D = init_D_from_txt('solway1.txt');

if take_map
    filename = sprintf('solway1_N=%d_alpha=%.4f_nsamples=%d_MAP.mat', N, h.alpha, nsamples);
else
    filename = sprintf('solway1_N=%d_alpha=%.4f_nsamples=%d_last.mat', N, h.alpha, nsamples);
end
disp(filename);

tic

for s = 1:N % for each simulated subject
    fprintf('infer H subject %d\n', s);

    if take_map
        [H, P] = sample_c(D, h, nsamples);
        [~,I] = max(P); % MAP H
        H = H(I);
    else
        [H, P] = sample_c(D, h, 1, nsamples);
    end
    chosen_H{s} = H;
end

toc;

save(filename, '-v7.3');

for s = 1:N % for each simulated subject
    fprintf('pick bus stop subject %d\n', s);

    H = chosen_H{s};

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


filename
save(filename);

%load('solway1.mat');

x = [-3 -1 -2 -3 -1 1 3 2 1 3];
y = [-3 -3 -2 -1 -1 1 1 2 3 3];
y = -y;

%{
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
%}

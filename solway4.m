% simulate experiment 4 from solway 2014

function filename = solway4(N, h, nsamples, take_map)

rng default;

sem = @(x) std(x) / sqrt(length(x));

if ~exist('N', 'var') || isempty(N)
    N = 35; % participants
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

D = init_D_from_txt('solway4.txt');

tasks = [
10 20;
5 22;
6 23;
8 18;
9 16;
19 12];

nexts = [
12 5;
6 10;
8 5;
6 9;
8 19;
20 9];

if take_map
    filename = sprintf('solway4_N=%d_alpha=%.4f_nsamples=%d_MAP.mat', N, h.alpha, nsamples);
else
    filename = sprintf('solway4_N=%d_alpha=%.4f_nsamples=%d_last.mat', N, h.alpha, nsamples);
end
disp(filename)

%load('solway4.mat');

tic

clear move;
for subj = 1:N % for each simulated subject
    fprintf('infer H subject %d\n', subj);

    if take_map
        [H, P] = sample_c(D, h, nsamples);
        [~,I] = max(P); % MAP H
        H = H(I);
    else
        [H, P] = sample_c(D, h, 1, nsamples);
    end
    chosen_H{subj} = H;
end

toc

save(filename, '-v7.3');

for subj = 1:N % for each simulated subject
    fprintf('HBFS subject %d\n', subj);

    H = chosen_H{subj};

    for t = 1:size(tasks,1)
        s = tasks(t,1);
        g = tasks(t,2);
        [path, hpath] = hbfs(s, g, H, D);
        move(subj, t) = path(2);

        % eps-greedy: choose random node w/ small prob 
        if rand() < 1 - h.eps
            move(subj, t) = datasample(find(D.G.E(s,:)), 1);
        end
    end
end

filename
save(filename);

%load('solway4.mat');
%load('solway4_1000.mat');

%{
c1 = 0;
c2 = 0;
n = 0;

synth = [];
for t = 1:size(tasks,1)
    c1 = c1 + sum(move(:,t) == nexts(t,1)); % count 1
    c2 = c2 + sum(move(:,t) == nexts(t,2)); % count 2
    synth = [synth ones(1,c1) zeros(1,c2)];
end
n = c1 + c2;
p = 2 * binocdf(min(c1,c2), n, 0.5);
y = binoinv([0.025 0.975], n, 0.5) / n;

figure;
m = c1/n;
se = sem(synth);
bar(m);
hold on;
errorbar(m, se);
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [y(1) y(1) y(2) y(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');
hold off;

fprintf('two-sided binomial test n = %d, p = %.4f\n', n, p);


%}

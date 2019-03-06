% simulate experiment 4 from solway 2014

function solway4(N, h, nsamples, take_map)

rng default;

sem = @(x) std(x) / sqrt(length(x));

if ~exist('N', 'var') || isempty(N)
    N = 35; % participants
end
if ~exist('h', 'var')
    h = init_hyperparams;
    h.alpha = 2;
end
if ~exist('nsamples', 'var')
    nsamples = 10000;
end
if ~exist('take_map', 'var')
    take_map = true;
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

%load('solway4.mat');


clear move;
for subj = 1:N % for each simulated subject
    fprintf('subject %d\n', subj);

    [H, P] = sample_c(D, h, nsamples);
    H_all{subj} = H;
    P_all{subj} = P;
    %H = H_all{subj};
    %P = P_all{subj};

    if take_map
        [~,I] = max(P); % MAP H
        H = H(I);
        map_H{subj} = H;
    else
        H = H(end); % last one
    end

    for t = 1:size(tasks,1)
        [path, hpath] = hbfs(tasks(t,1), tasks(t,2), H, D);
        move(subj, t) = path(2);
    end
end

filename = sprintf('solway4_N=%d_alpha=%.4f_nsamples=%d.mat', N, h.alpha, nsamples);
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

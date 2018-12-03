% simulate experiment 4 from solway 2014

clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 35; % participants
h.alpha = 5;


D = init_D_from_txt('solway4.txt');


for subj = 1:N % for each simulated subject
    fprintf('subject %d\n', subj);

    [H, P] = sample(D, h, 100);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);

    [path, hpath] = hbfs(10, 20, H, D);
    move(subj) = path(2);
end

save('solway3.mat');

load('solway3.mat');

c1 = sum(move ~= 5); % count 1
c2 = sum(move == 5); % count 2
n = c1 + c2;
p = 2 * binocdf(min(c1,c2), n, 1/3);
y = binoinv([0.025 0.975], n, 1/3);

figure;
m = c1/n;
se = sem(move ~= 5);
bar(m);
hold on;
errorbar(m, se);
line([0 2], [1/3 1/3], '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [y(1) y(1) y(2) y(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');
hold off;

fprintf('two-sided binomial test n = %d, p = %.4f\n', n, p);




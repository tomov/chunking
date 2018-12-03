% simulate experiment 4 from solway 2014

clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 35; % participants
h.alpha = 5;


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

load('solway4.mat');


clear move;
for subj = 1:N % for each simulated subject
    fprintf('subject %d\n', subj);

    %[H, P] = sample(D, h, 1000);
    %H_all{subj} = H;
    %P_all{subj} = P;
    H = H_all{subj};
    P = P_all{subj};

    [~,I] = max(P); % MAP H
    H = H(I);
    map_H{subj} = H;

    for t = 1:size(tasks,1)
        [path, hpath] = hbfs(tasks(t,1), tasks(t,2), H, D);
        move(subj, t) = path(2);
    end
end

save('solway4.mat');

load('solway4.mat');

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




% simulate RT experiment from Lynn et al. 2018
t
clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 78; % participants
h.alpha = 5;

nsteps = 1000;

D = init_D_from_txt('lynn.txt');

for i = 4:1 % important -- backwards
    A((D.G.E ^ i) > 0) = i; % how many steps to get to node X
end

for subj = 1:N % for each simulated subject
    fprintf('subject %d\n', subj);

    [H, P] = sample(D, h, 100);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);

    cross = [];
    % 2 = short violation
    % 3,4 = long violation
    viol = [ones(1,20)*2 ones(1,20)*3 ones(1,10)*4];
    viol = viol(randperm(length(viol));
    violations = ones(1, nsteps);
    violations(randsample(nsteps, 50)) = viol; % for each trial, how far to go (how many transitions to make)

    v = 1;
    path = [v];
    for i = 1:nsteps
        next = find(A(v,:) == violations(i));
        next = datasample(next, 1);
        cross = [cross H.c(v) == H.c(next)];
        v = next;
    end

    p_cross(subj, 1) = mean(cross(violations == 1));
    p_cross(subj, 2) = mean(cross(violations == 2));
    p_cross(subj, 3) = mean(cross(violations == 3 || violations == 4));
end


save('lynn.mat');

load('lynn.mat');


p_short = p_cross(:,2) - p_cross(:,1);
p_long = p_cross(:,3) - p_cross(:,1);

m1 = mean(p_short);
se1 = sem(p_short);
m2 = mean(p_long);
se2 = sem(p_long);

figure;
bar([m1 m2]);
hold on;
errorbar([m1 m2], [se1 se2]);
hold off;

[h, p, ci, stats] = ttest(p_short);
fprintf('short violations: t(%d) = %.2f, p = %.4f (one sample two-tailed t-test against 0)', stats.df, stats.tstat, p);

[h, p, ci, stats] = ttest(p_long);
fprintf('long violations: t(%d) = %.2f, p = %.4f (one sample two-tailed t-test against 0)', stats.df, stats.tstat, p);

[h, p, ci, stats] = ttest2(p_short, p_long);
fprintf('short vs. long violations: t(%d) = %.2f, p = %.4f (two sample two-tailed t-test against 0)', stats.df, stats.tstat, p);

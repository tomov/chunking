% simulate Experiment 6
%

clear all;
rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 95; % participants
h = init_hyperparams();
nsamples = 10000;
take_map = false;

if take_map
    filename = sprintf('mines10_alpha=%d_nsamples=%d_MAP.mat', h.alpha, nsamples);
else
    filename = sprintf('mines10_alpha=%d_nsamples=%d_last.mat', h.alpha, nsamples);
end

D = init_D_from_txt('mines10.txt');

tic

for s = 1:N % for each simulated subject
    fprintf('infer H for subject %d\n', s);

    if take_map
        [H, P] = sample_c(D, h, nsamples);
        [~,I] = max(P); % MAP H
        H = H(I);
    else
        [H, P] = sample_c(D, h, 1, nsamples);
    end
    chosen_H{s} = H;
end

toc

save(filename, '-v7.3');

choices = [];
for s = 1:N % for each simulated subject
    fprintf('HBFS for subject %d\n', s);

    H = chosen_H{s};

    [path, hpath] = hbfs(6, 1, H, D);

    if path(2) == 5
        choices = [choices 1];
        disp('yay!');
    else
        choices = [choices 0];
        disp('nay...');
    end
end


c1 = sum(choices);
c2 = sum(~choices);
m = mean(choices);
se = sem(choices);
n = N;
p = 2 * binocdf(min(c1, c2), n, 0.5);


fprintf('two-tailed binomial test c1 = %d (m = %.3f), n = %d, p = %.4f\n', c1, m, n, p);



filename
save(filename);

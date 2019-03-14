% simulate experiment 5
%

clear all;

sem = @(x) std(x) / sqrt(length(x));


N = 32; % participants
h = init_hyperparams();
nsamples = 10000;
take_map = false;

D = init_D_from_txt('mines.txt');

if take_map
    filename = sprintf('mines_alpha=%.4f_nsamples=%d_MAP.mat', h.alpha, nsamples);
else
    filename = sprintf('mines_alpha=%.4f_nsamples=%d_last.mat', h.alpha, nsamples);
end
filename

%{
tic

for s = 1:N % for each simulated subject
    fprintf('subject %d\n', s);

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
%}

load(filename);

h = init_hyperparams;

init_all_plots;

if take_map
    filename = sprintf('mines_alpha=%.4f_nsamples=%d_eps=%.4f_MAP.mat', h.alpha, nsamples, h.eps);
else
    filename = sprintf('mines_alpha=%.4f_nsamples=%d_eps=%.4f_last.mat', h.alpha, nsamples, h.eps);
end
filename


choices = [];
for s = 1:N % for each simulated subject
    fprintf('choice for subj %d', s);

    H = chosen_H{s};

    for i = 1:D.G.N
        pred(i) = H.theta(H.c(i));
    end

    if pred(3) > pred(7)
        choice = 1;
        disp('yay!');
    else
        choice = 0;
        disp('nay...');
    end

    % eps-greedy: choose randomly with small prob
    if rand() < 1 - h.eps
        choice = rand < 0.5;
    end

    choices = [choices choice];
end


c1 = sum(choices);
c2 = sum(~choices);
m = mean(choices);
se = sem(choices);
n = N;
p = 2 * binocdf(min(c1,c2), n, 0.5);


fprintf('two-tailed binomial test c1 = %d (m = %.3f), n = %d, p = %.4f\n', c1, m, n, p);



save(filename);

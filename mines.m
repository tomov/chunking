% simulate experiment 5
%

clear all;
rng default;

sem = @(x) std(x) / sqrt(length(x));


N = 32; % participants
h = init_hyperparams();
h.alpha = 1;
nsamples = 10000;
take_map = false;

D = init_D_from_txt('mines.txt');

if take_map
    filename = sprintf('mines_alpha=%d_nsamples=%d_MAP.mat', h.alpha, nsamples);
else
    filename = sprintf('mines_alpha=%d_nsamples=%d_last.mat', h.alpha, nsamples);
end
filename

choices = [];
for s = 1:N % for each simulated subject
    fprintf('subject %d\n', s);

    [H, P] = sample_c(D, h, nsamples);
    H_all{s} = H;
    P_all{s} = P;

    if take_map
        [~,I] = max(P); % MAP H
        H = H(I);
    else
        H = H(end);
    end
    map_H{s} = H;

    for i = 1:D.G.N
        pred(i) = H.theta(H.c(i));
    end

    if pred(3) > pred(7)
        choices = [choices 1];
        disp('yay!');
    else
        choices = [choices 0];
        disp('nay...');
    end
end


c = sum(choices);
m = mean(choices);
se = sem(choices);
n = N;
p = 1 - binocdf(c, n, 0.5);


fprintf('right-tailed binomial test m = %.3f, n = %d, p = %.4f\n', m, n, p);



save(filename);

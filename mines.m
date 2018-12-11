% simulate experiment 6
%

clear all;
rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 32; % participants
h = init_hyperparams();
nsamples = 1000;
D = init_D_from_txt('mines.txt');


choices = [];
for s = 1:N % for each simulated subject
    fprintf('subject %d\n', s);

    [H, P] = sample(D, h, nsamples);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);
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



filename = sprintf('mines_alpha=%d_nsamples=%d.mat', h.alpha, nsamples);
save(filename);

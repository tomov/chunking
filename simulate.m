

clear all;
rng default;

h.alpha = 1.5;
h.var_theta = 10;
h.theta_mean = 15;
h.var_mu = 10;
h.var_r = 5;

D = init_D_from_txt_dynamic('symmetric_rewards.txt');

%predict(D, h, M, burnin, lag, tau, ids)
tally = 0;
H_all = {};
post_all = {};
for i = 1:20
    [p, mu, H, post] = predict(D, h, 1000, 1, 1, 1, [2,7]); 
    %[p, mu, H, post] = predict(D, h, 2000, 1000, 10, 1, [2,7]); 
    H_all{i} = H;
    post_all{i} = post;
    if p(1) > p(2)
        tally = tally + 1;
    end
    disp(i);
end

% Get "START is not in the domain of the target or proposal distribution."
% for larger nsamples.
[~, I] = max(post);
figure;
plot_H(H(I), D);
disp(tally);

clear all;
rng default;

h.alpha = 1.5;
h.var_theta = 10;
h.theta_mean = 15;
h.var_mu = 10;
h.var_r = 5;

% Static rewards
D = init_D_from_txt_dynamic('static_rewards.txt');

n_subjects = 8;
tally = 0;
H_all = {};
post_all = {};
for i = 1:n_subjects
    % Get most likely H (based on posterior probabilities)
    [H, post] = sample(D, h, 500, 1, 1);
    [~, max_index] = max(post);
    H_max = H(max_index);
    % Use hierarchical BFS to predict the path taken given H
    [path, hpath] = hbfs(6, 1, H_max, D);
    % If the simulated test subject went right (i.e. to node 3),
    % increment the tally.
    if path(2) == 5
        tally = tally + 1;
    end
    % Store the hierarchies
    H_all{i} = H;
    post_all{i} = post;
    disp(i);
end
disp(tally)
[~, I] = max(post);
figure;
plot_H(H(I), D);
% Dynamic rewards
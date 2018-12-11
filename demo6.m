clear all;
rng default;

% Set hyperparameters
h.alpha = 1.5;
h.var_theta = 10;
h.theta_mean = 15;
h.var_mu = 10;
h.var_r = 5;

% Read data
D = init_D_from_txt_dynamic('static_rewards.txt');
n_subjects = 6;

% Dynamic rewards
n_trials = 100;
tally = 0;
H_all = {};
post_all = {};
H_prev = [];
for i = 1:n_subjects
    disp(i);
    for t = 1:n_trials
        % Redraw rewards
        if rand <= 0.2
            new_rewards = rand([1 3])*30;
            for j = 1:3
                D.r{j} = [new_rewards(1)];
            end
            for j = 4:6
                D.r{j} = [new_rewards(2)];
            end
            for j = 7:10
                D.r{j} = [new_rewards(3)];
            end
        end
        % Sample
        [H, post] = sample(D, H_prev, h, 10, 1, 1);
        % Store the hierarchies
        H_all{i, t} = H;
        post_all{i, t} = post;
        disp(t);
        H_prev = H(end);
    end
    % Get most likely H (based on posterior probabilities)
    [~, max_index] = max(post);
    H_max = H(max_index);
    % Use hierarchical BFS to predict the path taken given H
    [path, hpath] = hbfs(6, 1, H_max, D);
    % If the simulated test subject went right (i.e. to node 3),
    % increment the tally.
    if path(2) == 5
        tally = tally + 1;
    end
    H_prev = [];
    disp(tally/i);
end
disp(tally)
[~, I] = max(post);
figure;
plot_H(H(I), D);
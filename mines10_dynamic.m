% simulate experiment 6 but using dynamic rewards
%

clear all;
rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 95; % participants
h = init_hyperparams();
nsamples = 1000;
D = init_D_from_txt('mines10.txt');

n_trials = 100;


choices = [];
for s = 1:N % for each simulated subject
    fprintf('subject %d\n', s);

	H_prev = [];
    H = [];
    P = [];
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
        [H_cur, P_cur] = sample(D, h, nsamples / n_trials, 1, 1, H_prev);
        H = [H H_cur];
        P = [P P_cur];

        disp(t);
        % Store the hierarchies
        H_prev = H(end);
    end

    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);
    map_H{s} = H;

    [path, hpath] = hbfs(6, 1, H, D);

    if path(2) == 5
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



filename = sprintf('mines10_dynamic_alpha=%d_nsamples=%d.mat', h.alpha, nsamples);
save(filename);

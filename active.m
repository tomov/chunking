% experiment 7 simulation 
%


clear all;
rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 40; % participants
h = init_hyperparams();
h.alpha = 1;
burnin = 5000; % TODO x 10
lag = 100; % TODO x 10
nsamples = 50;

filename = sprintf('active_alpha=%d_nsamples=%d_burnin=%d_lag=%d.mat', h.alpha, nsamples);



D = init_D_from_txt('active_1.txt');

choices = [];
for subj = 1:N
    subj

    [H, P] = sample_c(D, h, nsamples, burnin, lag);
    H_all{subj} = H;

    entropy = [];
    for i = 1:length(D.G.hidden_edges)
        u = D.G.hidden_edges(i, 1);
        v = D.G.hidden_edges(i, 2);

        % D, (u,v) in E
        D_u_v = D;
        % unhide
        D_u_v.G.hidden_edges(i,:) = [];
        D_u_v.G.hidden_E(u,v) = 0;
        D_u_v.G.hidden_E(v,u) = 0;
        %_u_v set to 1
        D_u_v.G.E(u,v) = 1;
        D_u_v.G.E(v,u) = 1;

        % D, (u,v) not in E
        D_not_u_v = D;
        % unhide
        D_not_u_v.G.hidden_edges(i,:) = [];
        D_not_u_v.G.hidden_E(u,v) = 0;
        D_not_u_v.G.hidden_E(v,u) = 0;
        % u_v set to 0
        D_not_u_v.G.E(u,v) = 0;
        D_not_u_v.G.E(v,u) = 0;

        % H(H|D,a) = H(H|D,(u,v) in E) * Pr[(u,v) in E|D]
        %          + H(H|D,(u,v) not in E) * Pr[(u,v) not in E|D]
        entropy(subj,i) = approx_entropy(H, D_u_v, h) * Pr_edge(u, v, H, D, h) + ...
        approx_entropy(H, D_not_u_v, h) * Pr_not_edge(u, v, H, D, h);

        entropy(subj,:)
    end

    % TODO flip randomly when indifferent
    [~,a] = max(entropy(subj,:));
    choices(subj) = a;
end

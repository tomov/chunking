% experiment 7 simulation 
%


clear all;
rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 40; % participants
h = init_hyperparams();
h.alpha = 2;
burnin = 500; % TODO x 10
lag = 10; % TODO x 10
nsamples = 100;

resample = true;

filename = sprintf('active_alpha=%d_nsamples=%d_burnin=%d_lag=%d.mat', h.alpha, nsamples);

graph_files = {'active_1.txt', 'active_2.txt', ...
            'active_3.txt', 'active_4.txt', ...
            'active_5.txt', 'active_6.txt', ...
            'active_7.txt'};
exp_choices = [1 1 1 1 2 2 2];

graph_files =  {'active_6.txt'};
exp_choices = [2];


for g = 1:length(graph_files)

    D = init_D_from_txt(graph_files{g});
    D_all{g} = D;

    choices{g} = [];
    entropy = [];
    for subj = 1:N
        subj

        [H, logP] = sample_c(D, h, nsamples, burnin, lag);
        %H_all{g, subj} = H;

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

            % need bridges to compute this TODO sample
            for l = 1:length(H)
                H(l).b = get_H_b(H(l), D);
            end
            Pr_u_v = Pr_edge(u, v, H, D, h);
            Pr_not_u_v = Pr_not_edge(u, v, H, D, h);
            H = rmfield(H, 'b'); % b/c logpost_c doesn't like it 

            if resample
                % draw new samples after observation
                H_u_v = sample_c(D_u_v, h, nsamples, burnin, lag);
                H_not_u_v = sample_c(D_not_u_v, h, nsamples, burnin, lag);
            else
                % compute entropy using old samples
                H_u_v = H;
                H_not_u_v = H;
            end

            % H(H|D,a) = H(H|D,(u,v) in E) * Pr[(u,v) in E|D]
            %          + H(H|D,(u,v) not in E) * Pr[(u,v) not in E|D]
            entropy(subj,i) = approx_entropy(H, D_u_v, h) * Pr_u_v + ...
            approx_entropy(H, D_not_u_v, h) * Pr_not_u_v;

            %entropy(subj,:)
        end

        % TODO flip randomly when indifferent
        [~,a] = min(entropy(subj,:));
        choices{g}(subj) = a;
    end

    m(g) = mean(choices{g} == exp_choices(g));
    se(g) = sem(choices{g} == exp_choices(g));
end


figure;


for g = 1:length(graph_files)
    subplot(2, length(graph_files), g);
    plot_D(D_all{g});
end

subplot(2,1,2);
bar(m);
hold on;
errorbar(m, se, 'linestyle', 'none');
plot([0.5 4.5], [0.5 0.5], '--', 'color', [0.5 0.5 0.5]);
plot([4.5 7.5], [1/3 1/3], '--', 'color', [0.5 0.5 0.5]);

save(filename);

% simulate experiment 2 from solway 2014

clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 10; % participants
h.alpha = 5;

ntasks = 50; 
null_iters = 1000;

D = init_D_from_txt('solway2.txt');


for subj = 1:N % for each simulated subject
    fprintf('subject %d\n', subj);

    [H, P] = sample(D, h, 100);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);

    H = populate_H(H, D); % fill up bridges

    for i = 1:ntasks
        s = datasample(1:9, 1);
        g = datasample(11:19, 1);
        [path,~,~,b] = hbfs(s, g, H, D);
        b = b(2:end-1); % all bridge nodes on hierarchical path (exclude s and g)

        if isempty(b)
            b = path(2); % if no bridges on path, pick next node on path as "first that comes to mind"
        end

        spath = bfs(s, g, D.G.E); % find actual shortest path
        spath = spath(2:end-1);

        loc(subj,i) = b(1); % first bridge node on path
        corr(subj,i) = ismember(s, spath); % only count if on actual shortest path ("correct")
        
        for j = 1:null_iters
            null{j}(subj,i) = datasample(spath, 1);
        end
    end

    p(s) = mean(loc(s, corr(s,:)) == 10);
end


for j = 1:null_iters
    null_p(j) = mean(null{j}(corr(:)) == 10);
end
null_p = sort(null_p);
lcb = null_p(length(null_p) * 0.025);
ucb = null_p(length(null_p) * 0.975);

m = mean(p);
se = sem(p);

figure;
hold on;
bar(m);
hold on;
errorbar(m, se);
line([0 2], [mean(null_p) mean(null_p)], '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [lcb lcb ucb ucb], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');
hold off;

pos = find(m <= null_p);
pos = pos(end);
fprintf('MC test (%d samples from null), p = %.4f\n', null_iters, pos / length(null_p));

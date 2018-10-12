sem = @(x) std(x) / sqrt(length(x));

s = [];
g = [];
len = [];
group = [];

subj_group = [];
for subj = 1:size(data,1)
    phase = 2;
    for i = 1:length(data(subj, phase).s)
        s = [s; data(subj, phase).s(i)];
        g = [g; data(subj, phase).g(i)];
        len = [len; data(subj, phase).len(i)];
        group = [group; data(subj, phase).group(i)];
    end
    subj_group = [subj_group; data(subj,1).group(1)];
end


start = [2 4 6 7];
goal = [6 7 2 4];

figure;

for t = 1:length(start)
    A = group == 1 & s == start(t);
    B = group == 2 & s == start(t);

    m = [mean(len(A)), mean(len(B))];
    se = [sem(len(A)), sem(len(B))];

    [h, p, ci, stat] = ttest2(len(A), len(B));

    subplot(2,2,t);
    title(sprintf('%d -> %d: p = %.3f, t(%d) = %.3f', start(t), goal(t), p, stat.df, stat.tstat));
    hold on
    bar(1:2,m);
    errorbar(1:2,m,se,'.');
    xticklabels({'', 'A', 'B'});
    ylim([4 5]);

    if t == 1
        ylabel('action chunking')
    elseif t == 3
        ylabel('state chunking')
    end
end


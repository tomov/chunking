sem = @(x) std(x) / sqrt(length(x));

s = [];
g = [];
len = [];
group = [];
dir = [];

subj_group = [];
for subj = 1:size(data,1)
    phase = 2;
    for i = 1:length(data(subj, phase).s)
        s = [s; data(subj, phase).s(i)];
        g = [g; data(subj, phase).g(i)];
        len = [len; data(subj, phase).len(i)];
        dir = [dir; data(subj, phase).path{i}(2)];
        group = [group; data(subj, phase).group(i)];
    end
    subj_group = [subj_group; data(subj,1).group(1)];
end


start = [1 6 2 7];
goal = [6 1 7 2];

figure;

for t = 1:length(start)
    which = s == start(t);
    move = dir(which);
    m = unique(move); % kinds of second state (moves)
    assert(length(m) == 2);
    c1 = sum(move == m(1)); % count 1
    c2 = sum(move == m(2)); % count 2
    d = abs(c1 - c2);
    n = sum(which);
    p = 2 * binopdf((n - d) / 2, n, 0.5);

    subplot(2,2,t);
    bar(1:2, [c1 c2]);
    xticklabels({num2str(m(1)), num2str(m(2))});
    title(sprintf('%d -> %d: p = %.3f (d = %d, n = %d)', start(t), goal(t), p, d, n));
    %ylim([4 5]);

    if t == 1
   %     ylabel('action chunking')
    elseif t == 3
   %     ylabel('state chunking')
    end
end


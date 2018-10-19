% simulate data from chunking experiment

clear all;

%{
h.alpha = 1.5;

% load data
D = init_Ds_from_data('exp/results/subway10');
for i = 1:length(D)
    [samples, post] = sample(D(i), h, 1000);
    for j = 1:length(samples)
        H(i,j) = samples(j);
        P(i,j) = post(j);
    end
end

save demo5.mat;
%}

load demo5.mat;

%{
% hack sanity check -- make them all like D(1)
for i = 1:length(D)
    D(i) = D(1);
    H(i,:) = H(1,:);
    P(i,:) = P(1,:);
end
%}


% plot 5 per subject
%{
figure;
k = 5;
l = min(length(D), 10);
for i = 1:l
    post = P(i,:);
    [~,I] = maxk(post, k); % MAP k
    %I = length(post) - k : length(post); % last k 
    for j = 1:k
        subplot(l,k, (i-1)*k+j);
        plot_H(H(i,I(j)), D(i));
        if j == ceil(k/2);
            %ylabel(D(i).name);
            title(D(i).name);
        end
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
end
%}

% plot all subjects
figure;
rows = 5;
cols = 8;
s = 1;
for i = 1:rows
    for j = 1:cols
        if s > length(D)
            continue;
        end
        post = P(s,:);
        [~,I] = max(post); % MAP

        subplot(rows, cols, s);
        plot_H(H(s,I), D(s));
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        s = s + 1;
    end
end




% compute stats
% TODO dedupe w/ analyze_data.m

start = [6 7 1 2];
goal = [1 2 6 7];
nexts = [
5 7;
6 8;
2 10;
1 3
];

figure;


for t = 1:length(start)
    s = start(t);
    g = goal(t);
    clear move;
    for i = 1:length(D)
        [~,I] = maxk(P(i,:), 1);
        [path, hpath] = hbfs(s, g, H(i,I(1)), D(i));
        move(i) = path(2);
    end

    % from analyze_data.m
    m = nexts(t,:);
    c1 = sum(move == m(1)); % count 1
    c2 = sum(move == m(2)); % count 2
    d = abs(c1 - c2);
    n = length(D);
    p = 2 * binopdf((n - d) / 2, n, 0.5);

    subplot(2,2,t);
    bar(1:2, [c1 c2]);
    xticklabels({num2str(m(1)), num2str(m(2))});
    title(sprintf('%d -> %d: p = %.3f (d = %d, n = %d)', start(t), goal(t), p, d, n));
    %ylim([4 5]);

    if t == 1
        ylabel('state chunking')
    elseif t == 3
        ylabel('action chunking / S-A')
    end
end


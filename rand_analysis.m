% does the inferred clustering predict subject's test choices?

%analyze_data;

%load demo5.mat

assert(size(data,1) == length(D));

pred = [5 8 2 1];
c = [1 1 1 2 2 2 3 3 3 3];
rands = [];
sims = [];

for subj = 1:size(data,1)
    move = [dir(s == start(1) & ord == 1 & s_id == subj) ...
            dir(s == start(2) & ord == 1 & s_id == subj) ...
            dir(s == start(3) & ord == 1 & s_id == subj) ...
            dir(s == start(4) & ord == 1 & s_id == subj) ...
    ];
    move
    sim = sum(move == pred);
    sims = [sims; sim];

    [~, I] = max(P(subj,:));
    Hc = H(subj,I).c;
    randidx = RandIndex(c, Hc);
    rands = [rands; randidx];
end

[rho, p] = corr(rands, sims, 'type', 'Spearman');
rho
p

% simulate experiment 1 from solway 2013

clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 40; % participants
h.alpha = 5;

D = init_D_from_txt('solway1.txt');


for s = 1:N % for each simulated subject
    fprintf('subject %d\n', s);

    [H, P] = sample(D, h, 100);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);

    H = populate_H(H, D); % fill up bridges

    b = zeros(1,D.G.N);
    for i = 1:length(H.cnt)
        for j = i+1:length(H.cnt)
            bridge = H.b{i,j};
            if ~isempty(bridge)
                b(bridge(1)) = 1;
                b(bridge(2)) = 1;
            end
        end
    end

    b = find(b);
    if isempty(b)
        % only 1 cluster -> all nodes are fine
        b = 1:D.G.N;
    end
    loc(s,:) = datasample(b, 3); % pick 3 at random (with replacement)
end


save('solway1.mat');

load('solway1.mat');

H.c = [1 1 1 1 1 2 2 2 2 2];
h = plot_H(H, D);
for i = 1:D.G.N
    highlight(h, i, 'MarkerSize', sum(loc(:) == i));
end

function logp = loglik_task(H, D, h, s, g)
    % log P(task|G,H) 
    %

    logp = 0;

    logp = logp + log(1 / D.G.N);

    P = ones(1, D.G.N);
    P(H.c ~= H.c(s)) = H.tp;
    logp = logp + log(P(g)) - log(sum(P));

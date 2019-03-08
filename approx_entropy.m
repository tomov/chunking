function ent = approx_entropy(H, D, h)

    % H(H|D) = - sum P(H|D) log P(H|D)
    %     ~= - 1/M sum log p(m)
    %

    ent = 0;
    M = length(H);

    for i = 1:M
        logp(i) = logpost_c(H,D,h);
    end
    logp = lop - logsumexp(logp);

    for i = 1:M
        ent = ent + logp(i);
    end

    ent = - ent / M;

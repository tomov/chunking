function ent = approx_entropy(H, D, h)

    ent = 0;
    M = length(H);

    for i = 1:M
        logp(i) = logpost_c(H(i),D,h);
    end
    logp = logp - logsumexp(logp);
    p = exp(logp);
    p = p / sum(p); % just in case

    ent = - sum(p .* logp);

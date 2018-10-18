function logp = logprior(H, D, h)
    %
    % P(H)
    %

    % chunks
    logp = 0;
    cnt = zeros(1, max(H.c));
    cnt(H.c(1)) = 1;
    for i = 2:D.G.N
        c = H.c(i);
        if cnt(c) == 0
            logp = logp + log(h.alpha) - log(sum(cnt) + h.alpha);
        else
            logp = logp + log(cnt(c)) - log(sum(cnt) + h.alpha);
        end
        cnt(c) = cnt(c) + 1;
    end

    % p's
    logp = logp + log(betapdf(H.p,1,1)) + log(betapdf(H.q,1,1)) + log(betapdf(H.tp,1,1)) + log(betapdf(H.hp,1,1)); % TODO const
end


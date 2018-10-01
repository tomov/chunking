function logp = logprior(H, D, h)
    %
    % P(H)
    %

    logp = 0;
    cnt = zeros(size(H.cnt));
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
    assert(isequal(cnt, H.cnt));

    logp = logp + log(betapdf(H.p,1,1)) + log(betapdf(H.q,1,1)); % TODO const

    % TODO bridges

end


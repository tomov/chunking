% P(H)
%
function logp = logprior(H, D, h)
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

    for k = 1:H.N
        for l = 1:k-1
            if H.E(k,l)
                logp = logp + log(H.hp);
            else
                logp = logp + log(1 - H.hp);
            end
        end
    end
     
    % TODO bridges

    logp = logp + log(betapdf(H.p,1,1)) + log(betapdf(H.q,1,1) + log(betapdf(H.hp,1,1))) + log(betapdf(H.tp,1,1)); % TODO const
end


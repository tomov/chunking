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

    % for each cluster
    for k = 1:length(H.c)
        % account for impact of theta on posterior
        % below probability is Pr that the particular value of theta was
        % drawn given that it was drawn from a normal dist w mu = 0, var =
        % 100 ; Pr(theta_k = x | rest of H) = normpdf(x; 0, 100)
        logp = logp + log(normpdf(H.theta(k), h.theta_mean, h.var_theta));
    end
    
    for i = 1:D.G.N
        logp = logp + log(normpdf(H.mu(i), H.theta(H.c(i)), h.var_mu));
    end

    if isinf(logp)
        logp = -1e100;
    end
    
    % TODO bridges

end


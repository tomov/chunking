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
    %logp = logp + log(betapdf(H.p,1,1)) + log(betapdf(H.q,1,1)) + log(betapdf(H.tp,1,1)) + log(betapdf(H.hp,1,1)); % TODO const; also removed b/c cluster MATLAB sucks


    % cluster rewards
    for k = 1:length(H.c) % TODO bug? max(H.c)?
        % account for impact of theta on posterior
        % below probability is Pr that the particular value of theta was
        % drawn given that it was drawn from a normal dist w mu = 0, var =
        % 100 ; Pr(theta_k = x | rest of H) = normpdf(x; 0, 100)
        logp = logp + log(normpdf(H.theta(k), h.theta_mean, h.std_theta));
    end
   
    % state rewards 
    for i = 1:D.G.N
        logp = logp + log(normpdf(H.mu(i), H.theta(H.c(i)), h.std_mu));
    end

    % prevent -Infs = impossible events; equivalent to using a
    % Gaussian + uniform mixture
    if isinf(logp)
        logp = -1e100;
    end
end


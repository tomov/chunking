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
    assert(immse(get_H_cnt(H,D), cnt) < 1e-9);

    % p's
    %logp = logp + log(betapdf(H.p,1,1)) + log(betapdf(H.q,1,1)) + log(betapdf(H.tp,1,1)) + log(betapdf(H.hp,1,1)); % TODO const; also removed b/c cluster MATLAB sucks


    % TODO BUG these are fucked in matlab b/c of new clusters; comment out for now. they're fine in the C version
    % cluster rewards
    for k = 1:length(H.theta)
        % account for impact of theta on posterior
        % below probability is Pr that the particular value of theta was
        % drawn given that it was drawn from a normal dist w mu = 0, var =
        % 100 ; Pr(theta_k = x | rest of H) = normpdf(x; 0, 100)
        if k <= length(cnt) && cnt(k) > 0
            % skip empty clusters -- they don't count
            p = log(normpdf(H.theta(k), h.theta_mean, h.std_theta));
            if isinf(p)
                % prevent -Infs = impossible events; equivalent to using a
                % Gaussian + uniform mixture
                % do it in a "soft" way so MCMC can recover one by one
                logp = logp + 1e-100;
            else
                logp = logp + p;
            end
            %fprintf('N(theta[%d]): logp += %e = %e\n', k, p, logp);
        end
    end
   
    % state rewards 
    for i = 1:D.G.N
        p = log(normpdf(H.mu(i), H.theta(H.c(i)), h.std_mu));
        if isinf(p)
            % prevent -Infs = impossible events; equivalent to using a
            % Gaussian + uniform mixture
            % do it in a "soft" way so MCMC can recover one by one
            logp = logp + 1e-100;
        else
            logp = logp + p;
        end
        %fprintf('N(mu[%d], theta[%d]): logp += %e = %e\n', i, H.c(i), p, logp);
    end

    assert(~isinf(logp));
end


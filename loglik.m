function logp = loglik(H, D, h)
    % P(D|H) = P(G,tasks|H) = P(tasks|G,H) P(G|H) 
    %

    % connectivity
    logp = 0;
    for i = 1:D.G.N
        for j = 1:i-1
            if H.c(i) == H.c(j)
                if D.G.E(i,j)
                    logp = logp + log(H.p);
                else
                    logp = logp + log(1 - H.p);
                end
            else
                if D.G.E(i,j)
                    logp = logp + log(H.p * H.q);
                else
                    logp = logp + log(1 - H.p * H.q);
                end
            end

            % TODO bridges
        end
    end

    %fprintf('at 1 -> %.6f\n', logp);

    % find chunk transitive closures
    c = unique(H.c);
    A = D.G.E;
    for i = 1:length(c)
        n = H.c == c(i);
        A(n,n) = closure(D.G.E(n,n));
    end

    % penalize disconnected chunks
    for i = 1:D.G.N
        for j = 1:i-1
            if H.c(i) == H.c(j) && ~A(i, j)
                logp = logp - 100; % TODO const...
            end
        end
    end

    %fprintf('at 2 -> %.6f\n', logp);

    % (hierarchical) edges
    N = max(H.c);
    E = get_H_E(H, D); % TODO prior or likelihood? who knows..
    cnt = get_H_cnt(H, D);
    for k = 1:N
        if cnt(k) > 0
            for l = 1:k-1
                if cnt(l) > 0
                    if E(k,l)
                        logp = logp + log(H.hp);
                        %fprintf(' %.4f for edge (%d %d)\n', log(H.hp), k, l);
                    else
                        logp = logp + log(1 - H.hp);
                        %fprintf(' %.4f for edge -(%d %d)\n', log(1 - H.hp), k, l);
                    end
                end
            end
        end
    end

    %fprintf('at 3 -> %.6f\n', logp);

    % bridges
    for k = 1:N % TODO bug? max(H.c)?
        if cnt(k) > 0
            for l = 1:k-1
                if cnt(l) > 0 && E(k,l)
                    logp = logp + log(1) - log(cnt(k)) - log(cnt(l));
                    logp = logp - log(H.p * H.q); % b/c the bridge is always there, but we penalized / overcounted the corresponding edge when accounting for the connectivity of G
                end
            end
        end
    end

    %fprintf('at 4 -> %.6f\n', logp);

    % tasks
    for i = 1:length(D.tasks.s)
        s = D.tasks.s(i);
        logp = logp + log(1 / D.G.N);

        g = D.tasks.g(i);
        P = ones(1, D.G.N);
        P(H.c ~= H.c(s)) = H.tp;
        logp = logp + log(P(g)) - log(sum(P));
    end


    %fprintf('at 5 -> %.6f\n', logp);


    % rewards
    for i = 1:D.G.N
        for obs = 1:length(D.r{i})
            % Pr(r = x | rest of H)
            logp = logp + log(normpdf( D.r{i}(obs), H.mu(i), h.std_r ));
        end
    end

    %fprintf('at 6 -> %.6f\n', logp);
    
    if isinf(logp)
        logp = -1e100;
    end
end


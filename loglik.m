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

    % tasks
    for i = 1:length(D.tasks.s)
        s = D.tasks.s(i);
        logp = logp + log(1 / D.G.N);

        g = D.tasks.g(i);
        P = ones(1, D.G.N);
        P(H.c ~= H.c(s)) = H.tp;
        logp = logp + log(P(g)) - log(sum(P));
    end
end


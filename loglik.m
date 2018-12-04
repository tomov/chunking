function logp = loglik(H, D, h)
    % P(D|H) = P(G|H)  
    %

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
    
    for i = 1:D.G.N
        for obs = 1:length(D.r{i})
            % Pr(r = x | rest of H)
            %save loglik.mat;
            logp = logp + log(normpdf( D.r{i}(obs), H.mu(i), h.var_r ));
            %fprintf("loglik");
            %disp(logp);
        end
    end
    
    if isinf(logp)
        logp = -1e100;
    end
end

function A = closure(E)
    % compute the transitive closure of the graph G
    %
    A = (eye(size(E)) + E) ^ size(E,1) > 0;
end

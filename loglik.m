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
    
    
    for i = 1:D.G.N
        for obs = 1:length(D.r{i})
            % Pr(r = x | rest of H)
            logp = logp + log(normpdf( D.r{i}(obs), H.mu(i), h.var_r ));
            %fprintf("loglik");
            %disp(logp);
        end
    end
end


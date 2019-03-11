function pr = Pr_not_edge(u, v, H, D, h)
    %
    % Pr[(u,v) not in E|D]
    % = sum over H: Pr[(u,v) not in E|H] * P[H|D]
    % ~= 1/M sum: Pr[(u,v) not in E|H^(m)]

    M = length(H);
    pr = 0;

    for i = 1:M
        k = H(i).c(u);
        l = H(i).c(v);
        if k == l
            pr = pr + (1 - H(i).p);
        elseif isequal(H(i).b{k,l}, [u v]) || isequal(H(i).b{k,l}, [v u])
            % it's a bridge => edge has to be there
            pr = pr + 0;
        else
            pr = pr + (1 - H(i).p * H(i).q);
        end
    end

    pr = pr / M;


function pr = Pr_not_edge(u, v, H, D, h)
    %
    % Pr[(u,v) not in E|D]
    % = sum over H: Pr[(u,v) not in E|H] * P[H|D]
    % ~= 1/M sum: Pr[(u,v) not in E|H^(m)]

    M = length(H);
    pr = 0;

    for i = 1:M
        if H(i).c(u) == H(i).c(v)
            pr = pr + (1 - H(i).p);
        else
            pr = pr + (1 - H(i).p * H(i).q);
        end
    end

    pr = pr / M;


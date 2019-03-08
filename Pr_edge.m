function pr = Pr_edge(u, v, H, D, h)
    %
    % Pr[(u,v) in E|D]
    % = sum over H: Pr[(u,v) in E|H] * P[H|D]
    % ~= 1/M sum: Pr[(u,v) in E|H^(m)]

    M = length(H);
    pr = 0;

    for i = 1:M
        if H(i).c(u) == H(i).c(v)
            pr = pr + H(i).p;
        else
            pr = pr + H(i).p * H(i).q;
        end
    end

    pr = pr / M;


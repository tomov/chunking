function E = get_H_E(H, D)

    % compute edges of H (deterministically based on H.c and D.G.E)

    N = max(H.c);
    E = zeros(N);
    for i = 1:D.G.N
        for j = 1:i-1
            if H.c(i) ~= H.c(j) && D.G.E(i,j)
                fprintf('hierarchical edge (%d %d)\n', H.c(i), H.c(j));
                E(H.c(i), H.c(j)) = 1;
                E(H.c(j), H.c(i)) = 1;
            end
        end
    end

end

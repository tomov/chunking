function cnt = get_H_cnt(H, D)

    % get # of states in each chunk

    cnt = zeros(1, max(H.c));
    for i = 1:D.G.N
        c = H.c(i);
        cnt(c) = cnt(c) + 1;
    end

end

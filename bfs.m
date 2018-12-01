function path = bfs(s, g, E)

    N = size(E, 1);
    q = [s];
    used = zeros(1, N);
    prev = zeros(1, N);
    used(s) = 1;

    f = 1;
    while f <= length(q)
        i = q(f);
        for j = randperm(N) % go in random order
        %for j = 1:N
            if E(i,j) && ~used(j)
                q = [q j];
                used(j) = 1;
                prev(j) = i;
            end
        end
        f = f + 1;
    end

    path = [];
    i = g;
    while i > 0
        path = [i path];
        i = prev(i);
    end

end

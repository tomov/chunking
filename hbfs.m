function [path, hpath, paths] = hbfs(s, g, H, D)

    hpath = bfs(H.c(s), H.c(g), H.E);

    % get subgoal/bridge sequence b
    b = [s];
    for i = 1:length(hpath) - 1
        b = [b H.b{hpath(i), hpath(i+1)];
    end

    % low-level bfs for each subgoal
    b = [b g];
    path = [];
    paths = {};
    for i = 1:2:length(p) - 1
        assert(H.c(b(i)) == H.c(b(i+1)));
        E = D.G.E;
        E(H.c ~= H.c(b(i)), H.c ~= H.c(b(i))) = 0;
        p = bfs(b(i), b(i+1), E);
        paths = [paths; {p}];
        path = [path p];
    end

    % exclude repeats
    which = path(1:end-1) ~= path(2:end);
    path = path(logical([which 1])); 

end

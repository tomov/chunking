function paths = hamiltonians(D)
    paths = {};

    for s = 1:D.G.N
        path = [s];
        idx = [1];
        while ~isempty(path)

            cur = path(end);
            found = false;
            if length(path) == D.G.N
                paths = [paths; {path}];
                path
            else
                for i = idx(end):D.G.N
                    if D.G.E(cur,i) && ~ismember(i,path)
                        path = [path i];
                        idx(end) = i + 1;
                        idx = [idx 1];
                        found = true;
                        break;
                    end
                end
            end

            if ~found
                path = path(1:end-1);
                idx = idx(1:end-1);
            end
        end
    end

end

function get_adj(ix, x_offs, y_offs)

dr = [0 1 0 -1];
dc = [1 0 -1 0];

adj = [];

for r = 1:size(ix, 1)
    for c = 1:size(ix, 2)
        i = ix(r, c);
        a = [];
        for dir = 1:4
            nr = r + dr(dir);
            nc = c + dc(dir);
            if nr < 1 || nr > size(ix, 1) || nc < 1 || nc > size(ix, 2)
                a = [a -1];
            else
                a = [a ix(nr, nc)];
            end
        end
        a = [a r + y_offs c + x_offs];
        adj = [adj; a];

        a = [a(1) a(4) a(3) a(2) a(6) a(5)];
        fprintf('%d %d %d %d %d %d\n', a(1), a(2), a(3), a(4), a(5), a(6));
    end
end



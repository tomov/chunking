function b = get_H_b(H, D)

    % populate bridges (at random)
    
    N = max(H.c);
    cnt = get_H_cnt(H, D);
    E = get_H_E(H, D);

    % all potential bridges
    e = cell(N);
    for i = 1:D.G.N
        for j = 1:i-1
            if H.c(i) ~= H.c(j) && D.G.E(i,j)
                e{H.c(i), H.c(j)} = [e{H.c(i), H.c(j)}; i j];
                e{H.c(j), H.c(i)} = [e{H.c(j), H.c(i)}; j i];
            end
        end
    end

    % pick bridge at random
    b = cell(N);
    for k = 1:N
        if cnt(k) > 0
            for l = 1:k-1
                if cnt(l) > 0 && E(k,l)
                    assert(numel(e{k,l}) > 0);
                    i = randi(size(e{k,l}, 1));
                    b{k,l} = e{k,l}(i, :);
                    b{l,k} = e{l,k}(i, :);
                end
            end
        end
    end

end

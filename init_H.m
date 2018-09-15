function H = init_H(D, h)
    % Initializes a hierarchy H based on data D and hyperparameters h
    %
    % H.c = clusters / state chunks i.e. H.c(i) = state chunk of state i
    % H.cnt = number of states in each chunk
    % H.p =  TODO finish
    %

    % Chinese restaurant process
    H.c = [1];
    cnt = [1];
    for i=2:D.G.N
        c_new = find(mnrnd(1, [cnt h.alpha] / sum([cnt h.alpha])));
        H.c = [H.c c_new];
        if c_new > length(cnt)
            cnt = [cnt 1];
        else
            cnt(c_new) = cnt(c_new) + 1;
        end
    end
    H.cnt = cnt;

    H.p = betarnd(1,1); % TODO const 
    H.q = betarnd(1,1); % TODO const 
    H.tp = betarnd(1,1); % TODO const 

    H.N = length(cnt);
    H.hp = betarnd(1,1); % TODO const 
    H.E = zeros(H.N, H.N); % TODO sparse ?
    for k = 1:H.N
        for l = 1:k-1
            if rand < H.hp
                H.E(k,l) = 1;
                H.E(l,k) = 1;
            end
        end
    end

    % TODO bridges
end



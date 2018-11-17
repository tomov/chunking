function H = init_H(D, h)
    % Initializes a hierarchy H based on data D and hyperparameters h
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

    H.theta = normrnd(0, h.var_theta, [1, length(H.c)]); % sample each element of the vector from normal
    % initialize H.mu for mu's that we don't know
    H.mu = [];
    for i = 1:D.G.N
        H.mu(i) = normrnd(H.theta(H.c(i)), h.var_mu);
    end
   % H.r((1:D.G.N)) = normrnd(H.mu(H.c((1:D.G.N))), 0);
    H.p = betarnd(1,1); % TODO const 
    H.q = betarnd(1,1); % TODO const 
end



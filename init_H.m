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

    H.gammas = normrnd(len(H.c), 0, 100) % sample each element of the vector from normal
    % initialize H.mu for mu's that we don't know
    H.mu(H.which_mus_are_not_known) = normrnd(H.gammas(which_mus_are_not_known), 1000)
    H.p = betarnd(1,1); % TODO const 
    H.q = betarnd(1,1); % TODO const 
end



% proposal PMF for c_i
% inspired by Algorithm 5 from Neal 1998: MCMC for DP mixtures
%
% notice proposal distr q(x'|x) = q(x') i.e. it's independent of previous cluster assignment => implements Independent Metropolis-Hastings
function P = propP_c_i(i, H, D, h)
    cnt = get_H_cnt(H, D);
    cnt(H.c(i)) = cnt(H.c(i)) - 1;

    z = find(cnt == 0); % reuse empty bins -- notice this is legit b/c we're not reusing parameters; empty bin = new cluster; in fact, we have to take care of that in case c_i_old was the only one by itelf
    if isempty(z)
        cnt = [cnt h.alpha];
    else
        %cnt(z) = h.alpha; % TODO BUG WRONG -- this should be alpha / length(z) b/c we can have multiple empty bins; HOWEVER this actually ends up working in our favor b/c it is more lenient towards creating new clusters => keep it for now ; instead, use Algorithm 7 from Neal 1998
        cnt(z) = h.alpha / length(z); % notice all the empty bins have equal probability = alpha, but that's fine b/c it doesn't matter which one we use as the new cluster; we just have to make sure their total probability is not too high, otherwise effective alpha is greater 
    end
    P = cnt / sum(cnt);
end


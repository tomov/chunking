function [D, samples, post] = sample_graph_update(D, h, nwait_update, nsamples, nparticles, burnin, lag)
    if ~exist('nsamples', 'var')
        nsamples = 100;
    end

    if ~exist('burnin', 'var')
        burnin = 1; % no burn-in
    end

    if ~exist('lag', 'var')
        lag = 1;
    end

    if ~exist('nparticles', 'var')
        nparticles = 10;
    end
    
    if ~exist('nwait_update', 'var')
        nwait_update = 1;
    end
    
    % initialize
    W = zeros(nparticles, 1);
    
    if any(D.G.I(:))
        [H, W] = sample(D, h, nparticles, 1000, 10);
    else
        for i = 1:nparticles
            H(i) = init_H(D, h);
    %         [samples_t{i}, post_t{i}, H(i)] = sample_Hm(D, H(i), h, 1000, 1, 1);
            [~, ~, H(i)] = sample_Hm(D, H(i), h, 1000, 1, 1);
            W(i) = exp(loglik(H(i), D, h));
        end
    end

    W = W/sum(W);

    s =  size(D.updates);
    num_updates = s(1);
    for i = 1:num_updates
        new_edge = D.updates(i, :);
        D = update_D(D, new_edge);
        if mod(i, nwait_update) == 0
            for j = 1:length(H)
%                 [samples_t{j}, post_t{j}, H(j)] = sample_Hm(D, H(j), h, nsamples, 1, 10);
                [~, ~, H(j)] = sample_Hm(D, H(j), h, nsamples, 1, 10);
            end
        end
        W = update_weights(W, H, new_edge);
    end
    
    % sample from H
    pd = makedist('Multinomial','probabilities', W);
    r = random(pd);
%     samples = samples_t{r}; post = post_t{r};
    samples = H(r); post = 1;
end

function [samples, post, H] = sample_Hm(D, H, h, nsamples, burnin, lag)
    % Roberts & Rosenthal (2009)
    for n = 1:nsamples * lag + burnin
        
        % single MH step
        for i = 1:D.G.N
            logp = @(c_i) logpost_c_i(c_i, i, H, D, h);
            proprnd = @(c_i_old) proprnd_c_i(c_i_old, i, H, D, h);
            logprop = @(c_i_new, c_i_old) logprop_c_i(c_i_new, c_i_old, i, H, D, h);

            [c_i, accept] = mhsample(H.c(i), 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
            H = update_c_i(c_i, i, H);
        end

        logp = @(p) logpost_p(p, H, D, h);
        proprnd = @(p_old) proprnd_p(p_old, H, D, h);
        logprop = @(p_new, p_old) logprop_p(p_new, p_old, H, D, h);

        [p, accept] = mhsample(H.p, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop); % TODO adaptive
        H.p = p;

        [q, accept] = mhsample(H.q, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
        H.q = q;

        % TODO bridges
        samples(n) = H;
        post(n) = logpost(H,D,h);
    end

    post = post(burnin:lag:end);
    samples = samples(burnin:lag:end);
end

% P(H|D) up to proportionality constant
%
function logp = logpost(H, D, h)
    logp = loglik(H, D, h) + logprior(H, D, h);
end

% Update H.c(i) and counts
% TODO makes copy of H -- super slow...
%
function H = update_c_i(c_i, i, H)
    H.cnt(H.c(i)) = H.cnt(H.c(i)) - 1;
    H.c(i) = c_i;
    if c_i <= length(H.cnt)
        H.cnt(H.c(i)) = H.cnt(H.c(i)) + 1;
    else
        H.cnt = [H.cnt 1];
    end
end


% P(H|D) for updates of c_i
% i.e. with new c's up to c_i, the candidate c_i, then old c's after (and old rest of H)
%
function logp = logpost_c_i(c_i, i, H, D, h)
    H = update_c_i(c_i, i, H);
    logp = logpost(H, D, h);
end

% proposal PMF for c_i
% inspired by Algorithm 5 from Neal 1998: MCMC for DP mixtures
%
function P = propP_c_i(c_i_old, i, H, D, h)
    cnt = H.cnt;
    cnt(H.c(i)) = cnt(H.c(i)) - 1;
    z = find(cnt == 0); % reuse empty bins TODO is this legit?
    if isempty(z)
        cnt = [cnt h.alpha];
    else
        cnt(z) = h.alpha;
    end
    P = cnt / sum(cnt);
end

% propose c_i
%
function c_i_new = proprnd_c_i(c_i_old, i, H, D, h)
    P = propP_c_i(c_i_old, i, H, D, h);
    c_i_new = find(mnrnd(1, P));

    % TODO bridges
end

function [logP, P] = logprop_c_i(c_i_new, c_i_old, i, H, D, h) % TODO merge w/ proprnd
    P = propP_c_i(c_i_old, i, H, D, h);
    logP = log(P(c_i_new));
end


% P(H|D) for updates of p
%
function logp = logpost_p(p, H, D, h)
    H.p = p;
    logp = logpost(H, D, h);
end

% proposals for p; random walk 
%
function p_new = proprnd_p(p_old, H, D, h)
    while true % TODO can use universality of uniform inverse CDF thingy
        p_new = normrnd(p_old, 0.1); % TODO const TODO adaptive
        if p_new <= 1 && p_new >= 0
            break; % keep params within bounds
        end
    end
end

% account for truncating that keeps params within bounds 
%
function logp = logprop_p(p_new, p_old, H, D, h)
    Z = normcdf(1, p_old, 0.1) - normcdf(0, p_old, 0.1); % TODO consts TODO adaptive
    logp = log(normpdf(p_new, p_old, 1)) - log(Z);
end

function D = update_D(D, new_edge)
    i = new_edge(1); j = new_edge(2); exists = new_edge(3);

    D.G.E(i,j) = exists;
    D.G.E(j,i) = exists;
    D.G.I(i,j) = 1;
    D.G.I(j,i) = 1;
       
end

function W = update_weights(W, H, new_edge)
    assert(length(H) == length(W));
    for i = 1:length(W)
        W(i) = W(i) * exp(loglik_update(H(i), new_edge));
    end
    W = W/sum(W);
end

% assumes new_edge state is known
function logp = loglik_update(H, new_edge)
    i = new_edge(1); j = new_edge(2); exists = new_edge(3);
    if H.c(i) == H.c(j)
        if exists
            logp = log(H.p);
        else
            logp = log(1 - H.p);
        end
    else
        if exists
            logp = log(H.p * H.q);
        else
            logp = log(1 - H.p * H.q);
        end
    end
end

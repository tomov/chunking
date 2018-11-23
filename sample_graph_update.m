function [samples, post] = sample_graph_update(D, h, num_updates, V_updates, E_updates, nsamples, nparticles, burnin, lag)
    if ~exist('nsamples', 'var')
        nsamples = 10000;
    end

    if ~exist('burnin', 'var')
        burnin = 1; % no burn-in
    end

    if ~exist('lag', 'var')
        lag = 1;
    end

    if ~exist('nparticles', 'var')
        nparticles = 1;
    end
    
    nparticles = 1;
    num_updates = 1;

    for i = 1:nparticles
        H(i) = init_H(D, h);
    end
    
    W = ones(nparticles, 1)*1/nparticles;

    for i = 1:num_updates
        assert(V_updates(i) + length(E_updates{i}) > 0);
        D = update_D(D, V_updates(i), E_updates{i});
%       TODO
        if V_updates(i) > 0
            % assign clusters to new vertices
            H = update_H(D, H, num_vertices, h);
        end
        % sample
        for j = 1:length(H)
            [samples, post, H(j)] = sample_Hm(D, H(j), h, nsamples, burnin, lag);
%             W = update_weights(W);
        end

    end
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

function D = update_D(D, num_vertices, new_edges)
    % add vertices 
    adds = num_vertices;
    D.G.E = [D.G.E; zeros(adds, D.G.N)]
    D.G.N = D.G.N + adds;
    D.G.E = [D.G.E zeros(D.G.N, adds)]
    
    % add lines
    adds = length(new_edges);
    for k = 1:adds
        edge = new_edges{k};
        i = edge(1); j = edge(2);
        assert(D.G.N >= i);
        assert(D.G.N >= j);
        assert(D.G.E(i,j) == 0);
        assert(D.G.E(j,i) == 0);
        D.G.E(i,j) = 1;
        D.G.E(j,i) = 1;
    end
       
end

%TODO
function H = update_H(D, H, num_vertices, h)
    % Chinese restaurant process
    for m = 1:length(H)
        % assign clusters to new vertices
        cnt = H(m).cnt;
        for i = D.G.N-num_vertices:D.G.N
            c_new = find(mnrnd(1, [cnt h.alpha] / sum([cnt h.alpha])));
            H(m).c = [H(m).c c_new];
            if c_new > length(cnt)
                cnt = [cnt 1];
            else
                cnt(c_new) = cnt(c_new) + 1;
            end
        end
        H(m).cnt = cnt;
    end
end

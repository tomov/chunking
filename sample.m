function [samples, post] = sample(D, h, nsamples, burnin, lag, H)
    %
    % Draw samples from posterior P(H|D) using Metropolis-Hastings-within-Gibbs sampling.
    % hierarchy H = (c, p, q, p', p", E', V')
    % data D = (G, tasks)
    % graph G = (E, V)
    % tasks = (task_1, task_2 ...)
    % task = (s, g)
    %
    % Generative model:
    %
    % P(H):
    % state chunks = c ~ CRP
    % within-cluster density = p ~ Beta
    % across-cluster density = pq, q ~ Beta
    % H graph density = p' = hp ~ Beta
    % probability goal state is in different chunk from starting state = p" = tp ~ Beta
    % 
    % P(G|H):
    % E(i,j) ~ Bern(p) if c(i) == c(j)
    % E(i,j) ~ Bern(pq) if c(i) != c(j)
    % 
    % P(tasks|G,H) = product P(task|G,H)
    % P(task|G,H):
    % starting state = s ~ Cat(all vertices in G)
    % goal state = g ~ Cat(1,1,1,1... for all i s.t. c(i) == c(s), ... p", p", p"... for all i s.t. c(i) != c(s)
    %

    if ~exist('nsamples', 'var')
        nsamples = 10000;
    end

    if ~exist('burnin', 'var')
        burnin = 1; % no burn-in
    end

    if ~exist('lag', 'var')
        lag = 1;
    end

    if ~exist('H', 'var') || isempty(H)
        H = init_H(D, h);
    end

    % Roberts & Rosenthal (2009)
    for n = 1:nsamples * lag + burnin
        for i = 1:D.G.N
            logp = @(c_i) logpost_c_i(c_i, i, H, D, h);
            proprnd = @(c_i_old) proprnd_c_i(c_i_old, i, H, D, h);
            logprop = @(c_i_new, c_i_old) logprop_c_i(c_i_new, c_i_old, i, H, D, h);

            [c_i, accept] = mhsample(H.c(i), 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
            H.c(i) = c_i;
        end

        logp = @(p) logpost_p(p, H, D, h);
        proprnd = @(p_old) proprnd_p(p_old, H, D, h);
        logprop = @(p_new, p_old) logprop_p(p_new, p_old, H, D, h);

        [p, accept] = mhsample(H.p, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop); % TODO adaptive
        H.p = p;

        logp = @(q) logpost_q(q, H, D, h);
        [q, accept] = mhsample(H.q, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
        H.q = q;

        logp = @(tp) logpost_tp(tp, H, D, h);
        [tp, accept] = mhsample(H.tp, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
        H.tp = tp;

        logp = @(hp) logpost_hp(hp, H, D, h);
        [hp, accept] = mhsample(H.hp, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
        H.hp = hp;

		logprop_thetamu = @(p_new, p_old) logprop_unbounded(p_new, p_old, H, D, h);

		% thetas
        for k = 1:length(H.c)
            proprnd_thetak = @(p_old) proprnd_unbounded(p_old, H, D, h);
            logtheta = @(theta_k) logpost_theta(theta_k, k, H, D, h);
            [theta_k, accept] = mhsample(H.theta(k), 1, 'logpdf', logtheta, 'proprnd', proprnd_thetak, 'logproppdf', logprop_thetamu);
            H.theta(k) = theta_k;
        end
        % do sampling for mu's
        for i = 1:D.G.N
            proprnd_mui = @(p_old) proprnd_unbounded(p_old, H, D, h);
            logmu = @(mu_i) logpost_mu(mu_i, i, H, D, h);
            %save debug.mat;
            [mu_i, accept] = mhsample(H.mu(i), 1, 'logpdf', logmu, 'proprnd', proprnd_mui, 'logproppdf', logprop_thetamu);
            H.mu(i) = mu_i; 
        end

        % TODO bridges

        samples(n) = H;
        post(n) = logpost(H,D,h);
    end

    samples = samples(burnin:lag:end);
end



% P(H|D) for updates of c_i
% i.e. with new c's up to c_i, the candidate c_i, then old c's after (and old rest of H)
%
function logp = logpost_c_i(c_i, i, H, D, h)
    H.c(i) = c_i;
    logp = logpost(H, D, h);
end

% proposal PMF for c_i
% inspired by Algorithm 5 from Neal 1998: MCMC for DP mixtures
%
function P = propP_c_i(c_i_old, i, H, D, h)
    cnt = get_H_cnt(H, D);
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

% P(H|D) for updates of q
%
function logp = logpost_q(q, H, D, h)
    H.q = q;
    logp = logpost(H, D, h);
end

% P(H|D) for updates of tp
%
function logp = logpost_tp(tp, H, D, h)
    H.tp = tp;
    logp = logpost(H, D, h);
end

% P(H|D) for updates of hp
%
function logp = logpost_hp(hp, H, D, h)
    H.hp = hp;
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
    logp = log(normpdf(p_new, p_old, 0.1)) - log(Z);
end


    
% P(H|D) for updates of theta
%
function logtheta = logpost_theta(theta_k, k, H, D, h)
    H.theta(k) = theta_k;
    logtheta = logpost(H, D, h);
end

% P(H|D) for updates of mu
%
function logmu = logpost_mu(mu_i, i, H, D, h)
    H.mu(i) = mu_i;
    logmu = logpost(H, D, h);
end


% proposals for mu, theta; random walk 
%
function p_new = proprnd_unbounded(p_old, H, D, h)
   % p_new = normrnd(p_old, exp(ls_s_k)); % TODO adaptive 
    p_new = normrnd(p_old, 1); % TODO const
end

% unbounded proposal f'n
% 
function logp = logprop_unbounded(p_new, p_old, H, D, h)
    logp = log(normpdf(p_new, p_old, 1));
end

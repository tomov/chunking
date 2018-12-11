function [samples, post] = sample(D, H, h, nsamples, burnin, lag)
    %
    % Draw samples from posterior P(H|D) using Metropolis-Hastings-within-Gibbs sampling.
    % hierarchy H = (c, p, q, p', p", E', V')
    % data D = (G)
    % graph G = (E, V)
    %
    % Generative model:
    %
    % P(H):
    % state chunks = c ~ CRP
    % within-cluster density = p ~ Beta
    % across-cluster density = pq, q ~ Beta
    % 
    % P(G|H):
    % E(i,j) ~ Bern(p) if c(i) == c(j)
    % E(i,j) ~ Bern(pq) if c(i) != c(j)
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

    if isempty(H)
        H = init_H(D, h);
    end
    % Adaptive sampling: l_i's (log of std dev for each latent variable)
    ls_mu = zeros(1, D.G.N);
    ls_theta = zeros(1, length(H.c));
    
    % Roberts & Rosenthal (2009)
    for n = 1:nsamples * lag + burnin
        
        %save wtf.mat
        
        for i = 1:D.G.N
            logp = @(c_i) logpost_c_i(c_i, i, H, D, h);
            proprnd = @(c_i_old) proprnd_c_i(c_i_old, i, H, D, h);
            logprop = @(c_i_new, c_i_old) logprop_c_i(c_i_new, c_i_old, i, H, D, h);
            %save cs.mat
            [c_i, accept] = mhsample(H.c(i), 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop);
            H = update_c_i(c_i, i, H);
        end

        logp = @(p) logpost_p(p, H, D, h);
        logq = @(q) logpost_q(q, H, D, h);
        
        proprnd = @(p_old) proprnd_p(p_old, H, D, h);
        logprop = @(p_new, p_old) logprop_p(p_new, p_old, H, D, h);
        
        logprop_thetamu = @(p_new, p_old) logprop_unbounded(p_new, p_old, H, D, h);

        [p, accept] = mhsample(H.p, 1, 'logpdf', logp, 'proprnd', proprnd, 'logproppdf', logprop); % TODO adaptive
        H.p = p;

        [q, accept] = mhsample(H.q, 1, 'logpdf', logq, 'proprnd', proprnd, 'logproppdf', logprop);
        H.q = q;
        
        % do sampling for thetas
        % for proposal function, just use logprop_p (Gaussian random walk;
        % maybe want a different random walk that takes larger jumps on 
        % each step, e.g. more like size = 1 instead of < 1).
        % Use same target distribution (still proportional to the same
        % posterior)
        for k = 1:length(H.c)
            proprnd_thetak = @(p_old) proprnd_unbounded(p_old, H, D, h, ls_theta(1, k));
            logtheta = @(theta_k) logpost_theta(theta_k, k, H, D, h);
            [theta_k, accept] = mhsample(H.theta(k), 1, 'logpdf', logtheta, 'proprnd', proprnd_thetak, 'logproppdf', logprop_thetamu);
            H.theta(k) = theta_k;
            % Update l_i's to maintain acceptance rate of 0.44
            d = min(0.01, sqrt(n));
            if accept > 0.44
                ls_theta(1,k) = ls_theta(1,k) + d;
            else
                ls_theta(1,k) = ls_theta(1,k) - d;
            end
        end
        % do sampling for mu's
        for i = 1:D.G.N
            proprnd_mui = @(p_old) proprnd_unbounded(p_old, H, D, h, ls_mu(1, k));
            logmu = @(mu_i) logpost_mu(mu_i, i, H, D, h);
            %save debug.mat;
            [mu_i, accept] = mhsample(H.mu(i), 1, 'logpdf', logmu, 'proprnd', proprnd_mui, 'logproppdf', logprop_thetamu);
            H.mu(i) = mu_i; 
            % Update l_i's to maintain acceptance rate of 0.44
            d = min(0.01, sqrt(n));
            if accept > 0.44
                ls_mu(1,i) = ls_mu(1,i) + d;
            else
                ls_mu(1,i) = ls_mu(1,i) - d;
            end
        end
        
        % TODO bridges

        samples(n) = H;
        post(n) = logpost(H,D,h);
    end

    samples = samples(burnin:lag:end);
    post = post(burnin:lag:end);
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

% P(H|D) for updates of q
%
function logq = logpost_q(q, H, D, h)
    H.q = q;
    logq = logpost(H, D, h);
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

% proposals for mu, theta; random walk 
%
function p_new = proprnd_unbounded(p_old, H, D, h, ls_s_k)
    p_new = normrnd(p_old, exp(ls_s_k)); % TODO const
end

% account for truncating that keeps params within bounds 
%
function logp = logprop_p(p_new, p_old, H, D, h)
    Z = normcdf(1, p_old, 0.1) - normcdf(0, p_old, 0.1); % TODO consts TODO adaptive
    logp = log(normpdf(p_new, p_old, 0.1)) - log(Z);
    %save logprop.mat
end

% 
function logp = logprop_unbounded(p_new, p_old, H, D, h)
    logp = log(normpdf(p_new, p_old, 1));
    %save logprop.mat
end


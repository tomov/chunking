% P(H|D) up to proportionality constant
%
function logp = logpost(H, D, h)
    logp = loglik(H, D, h) + logprior(H, D, h);
end


% likelihood of single observation at time t
% used for particle filtering

function [lik, out] = pf_lik(t, particle, D, h)

    out = -1; % dummy

    if t == 1
        % edge case: initialize with likelihood of graph and tasks, P(D|H)
        D.tasks.s = D.tasks.s(1);
        D.tasks.g = D.tasks.g(1);
        
        lik = (loglik_c(particle.H, D, h));
    else
        % just take task into account, P(task|D)
        s = D.tasks.s(t);
        g = D.tasks.g(t);

        lik = (loglik_task(particle.H, D, h, s, g));
    end



% see MCMC_*_choice in https://github.com/tomov/Context-Learning-Data-Analysis/tree/isl

function [lik, out] = MH_choice(t, particle, D, h)

    s = D.tasks.s(t);
    g = D.tasks.g(t);

    [path, hpath] = hbfs(s, g, particle.H, D);
    out = path(2);

    a = D.path{t}(2);

    if a == out
        lik = 1 - h.eps;
    else
        lik = h.eps; % TODO assumes only 2 actions available at each time
    end

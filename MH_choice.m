% see MCMC_*_choice in https://github.com/tomov/Context-Learning-Data-Analysis/tree/isl

function [lik, out] = MH_choice(t, particle, D, h)

    s = D.tasks.s(t);
    g = D.tasks.g(t);

    [path, hpath] = hbfs(s, g, particle.H, D);
    if length(path) < 2
        % happens for shitty (disconnected) hierarchies
        warning('path length < 2!')
        %particle.H
        out = -1;
    else
        out = path(2);
    end

    a = D.path{t}(2);

    if a == out
        lik = h.eps; % note eps means opposite of convention
    else
        lik = 1 - h.eps; % TODO assumes only 2 actions available at each time
    end

    %fprintf('%d vs %d; %.4f\n', a, out, lik);


% see MCMC_*_update in  https://github.com/tomov/Context-Learning-Data-Analysis/tree/isl

function particle = MH_update(t, particle, D, h, nsamples, T)

    % b/c D is full data
    D.tasks.s = D.tasks.s(1:t);
    D.tasks.g = D.tasks.g(1:t);

    particle.H = sample_c(D, h, 1, round(nsamples / T), 1, particle.H);

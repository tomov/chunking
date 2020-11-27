% see MCMC_*_init in https://github.com/tomov/Context-Learning-Data-Analysis/tree/isl

function particle = MH_init(D, h)

    % D: subject data, as from init_D_from_data()
    % h: init_hyperparams

    particle.H = sample_c(D, h, 1); % take 1 sample to start with

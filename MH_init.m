% see MCMC_*_init in https://github.com/tomov/Context-Learning-Data-Analysis/tree/isl

function particle = MH_init(D, h)

    % D: subject data, as from init_D_from_data()
    % h: init_hyperparams

    % b/c D is full data
    D.tasks.s = [];
    D.tasks.g = [];

    particle.H = sample_c(D, h, 1, 100, 1); % 100 MCMC iterations to find H that is not disconnected TODO should be kosher b/c like rejuvination? keeps posterior invariant?

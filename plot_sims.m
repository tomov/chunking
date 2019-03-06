% plot all simulations
close all;
clear all;

nsamples = 10000;
h = init_hyperparams();
h.alpha = 1;

filename{1} = schapiro([], h, nsamples, false);
fig_schapiro(filename{1});

filename{2} = solway1([], h, nsamples, false);
fig_solway1(filename{2});

filename{3} = solway2([], h, nsamples, false);
fig_solway2(filename{3});

filename{4} = solway4([], h, nsamples, false);
fig_solway4(filename{4});

filename{5} = lynn([], h, nsamples, false);
fig_lynn(filename{5});

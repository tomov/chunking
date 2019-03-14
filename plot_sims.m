% plot all simulations
close all;
clear all;

nsamples = 10000;
take_map = false;
h = init_hyperparams();

filename{1} = schapiro([], h, nsamples, take_map);
fig_schapiro(filename{1});

filename{2} = solway1([], h, nsamples, take_map);
fig_solway1(filename{2});

filename{3} = solway2([], h, nsamples, take_map);
fig_solway2(filename{3});

filename{4} = solway4([], h, nsamples, take_map);
fig_solway4(filename{4});

filename{5} = lynn([], h, nsamples, take_map);
fig_lynn(filename{5});

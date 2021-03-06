% explore grid of parameters -- lots of samples, 1 markov chain

N = 1;
h = init_hyperparams;
nsamples = 100000;

for alpha = 1 : 0.5 : 20
    h.alpha = alpha;
    disp(h)
    tic
    schapiro(N, h, nsamples);
    solway1(N, h, nsamples);
    solway2(N, h, nsamples);
    solway4(N, h, nsamples);
    lynn(N, h, nsamples);
    toc
end

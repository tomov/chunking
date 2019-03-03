% explore grid of parameters -- 10k samples, N chains

N = [];
h = init_hyperparams;
nsamples = 10000;

for alpha = 3 : 0.5 : 20
    h.alpha = alpha;
    disp(h)
    tic
    schapiro(N, h, nsamples, false);
    solway1(N, h, nsamples, false);
    solway2(N, h, nsamples, false);
    solway4(N, h, nsamples, false);
    lynn(N, h, nsamples, false);
    toc
end

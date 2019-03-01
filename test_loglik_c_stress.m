% stresstest loglik_c vs. loglik

files = {'schapiro_N=1_alpha=1_nsamples=100000.mat', ...
         'solway1_N=1_alpha=1.0000_nsamples=100000.mat', ...
         'solway2_N=1_alpha=1.0000_nsamples=100000.mat'};

for f = 1:length(files)
    file = files{f};
    load(file);

    file
    for i = 1:1000 % length(H_all{1})
        l1 = loglik(H_all{1}(i), D, h);
        l2 = loglik_c(H_all{1}(i), D, h);

        assert(abs(l1 - l2) < 1e-9);
    end
end

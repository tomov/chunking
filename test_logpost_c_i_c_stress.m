% stresstest proposal distribution 

files = {'schapiro_N=32_alpha=2.0000_nsamples=1000.mat',  ...
         'schapiro_N=1_alpha=1_nsamples=100000.mat', ...
         'solway1_N=1_alpha=1.0000_nsamples=100000.mat', ...
         'solway2_N=1_alpha=1.0000_nsamples=100000.mat'};

for f = 1:length(files)
    file = files{f};
    load(file);

    file
    for t = 1:1000 % length(H_all{1})
        H = H_all{1}(t);
        for i = 1:length(H_all{1}(t).c)
            P1 = propP_c_i(i, H, D, h);

            for c_i_new = 1:length(P1)
                l1 = logpost_c_i(c_i_new, i, H, D, h);
                l2 = logpost_c_i_c(c_i_new, i, H, D, h);

                assert(abs(l1 - l2) < 1e-9);
            end
        end
    end
end

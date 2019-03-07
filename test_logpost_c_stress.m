% stresstest logpost_c vs. logpost

%{
files = {'schapiro_N=32_alpha=2.0000_nsamples=1000.mat',  ...
         'schapiro_N=1_alpha=1_nsamples=100000.mat', ...
         'solway1_N=1_alpha=1.0000_nsamples=100000.mat', ...
         'solway2_N=1_alpha=1.0000_nsamples=100000.mat'};

for f = 1:length(files)
    file = files{f};
    load(file);

    file
    for i = 1:1000 % length(H_all{1})
        H = H_all{1}(i);
        l1 = logpost(H, D, h);
        l2 = logpost_c(H, D, h);

        assert(abs(l1 - l2) < 1e-9);
    end
end
%}


% test rewards

%load model_exp_6_samples=1000_alpha=1.0000_MAP.mat;

%H = pl(5).H{1};
%D = pl(5).D{1};
D = D_all(1);

for i = 1:length(H_all(1,:))
    H = H_all(1,i);
    l1 = logpost(H, D, h);
    l2 = logpost_c(H, D, h);

    assert(abs(l1 - l2) < 1e-9);
end

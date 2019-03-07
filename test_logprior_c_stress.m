% stresstest logprior_c vs. logprior

%files = {'schapiro_N=32_alpha=2.0000_nsamples=1000.mat',  ...
%         'schapiro_N=1_alpha=1_nsamples=100000.mat', ...
%         'solway1_N=1_alpha=1.0000_nsamples=100000.mat', ...
%         'solway2_N=1_alpha=1.0000_nsamples=100000.mat'};
%
%for f = 1:length(files)
%    file = files{f};
%    load(file);
%
%    file
%    for i = 1:1000 % length(H_all{1})
%        l1 = logprior(H_all{1}(i), D, h);
%        l2 = logprior_c(H_all{1}(i), D, h);
%
%        assert(abs(l1 - l2) < 1e-9);
%    end
%end


for i = 1:length(H_all{1})
    H = H_all{1}(i);
    l1 = logprior_agni(H, D, h);

    H_shit.c = H.c;
    H_shit.p = H.p;
    H_shit.q = H.q;
    H_shit.tp = H.tp;
    H_shit.hp = H.hp;
    H_shit.theta = H.theta;
    H_shit.mu = H.mu;
    l2 = logprior_c(H_shit, D, h);

    assert(abs(l1 - l2) < 1e-9);
end

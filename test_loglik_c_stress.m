% stresstest loglik_c vs. loglik

%files = {'schapiro_N=32_alpha=2.0000_nsamples=1000.mat',  ...
%         'schapiro_N=1_alpha=1_nsamples=100000.mat', ...
%         'solway1_N=1_alpha=1.0000_nsamples=100000.mat', ...
%         'solway2_N=1_alpha=1.0000_nsamples=100000.mat'};

%for f = 1:length(files)
%    file = files{f};
%    load(file);
%
%    file
%    for i = 1:1000 % length(H_all{1})
%        l1 = loglik(H_all{1}(i), D, h);
%        l2 = loglik_c(H_all{1}(i), D, h);
%
%        assert(abs(l1 - l2) < 1e-9);
%    end
%end

% test tasks
%{
load model_all_data_samples=40_MAP_alpha=2.0000.mat;

H = pl(2).H{4}(1,:);
D = pl(2).D{4}(1);

for i = 1:length(H)
    l1 = loglik(H(i), D, h);
    l2 = loglik_c(H(i), D, h);

    assert(abs(l1 - l2) < 1e-9);
end
%}

for i = 1:length(H_all{1})
    H = H_all{1}(i);
    l1 = loglik_agni(H, D, h);

    H_shit.c = H.c;
    H_shit.p = H.p;
    H_shit.q = H.q;
    H_shit.tp = H.tp;
    H_shit.hp = H.hp;
    H_shit.theta = H.theta;
    H_shit.mu = H.mu;
    l2 = loglik_c(H_shit, D, h);

    assert(abs(l1 - l2) < 1e-9);
end

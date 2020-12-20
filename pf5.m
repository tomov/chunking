% this is the sameas particle filter 4,exceptcopy from particle folder 2
% simple particle filter for experiment 3, with rejuvination
% note this is NOT ISL, but approximately ideal observer
% c/p pf1.m
%
% for each particle,
%   H ~ P(H)
%   for each trial,
%     w = P(D|H)   % i.e. likelihood weighing, with full data D (so far)
%     H ~ MCMC(H)
%

clear all;
close all;

h = init_hyperparams;
nsamples = 1000;
num_particles = 1000;

index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..; ORDER CRUCIAL

sem = @(x) std(x) / sqrt(length(x));


% from model_all_data

filename = sprintf('mat/pf5_alpha=%.4f_nsamples=%d_div_eps=%.4f_last_np=%d.mat', h.alpha, nsamples, h.eps, num_particles);
[~, name] = system('hostname');
if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
    % local
    filename
else
    % cannon
    filename = fullfile(getenv('MY_SCRATCH'), 'chunking', filename);
    filename
end


D = init_Ds_from_data('exp/results/exp_v2_3_subway10_unlearn_circ', true);
D_full = D;


for subj = 1:length(D) % for each subject
    fprintf('infer H: subject %d\n', subj);

    tic

    T = length(D_full(subj).tasks.s);
    assert(T == 246);
    assert(T == length(D(subj).tasks.g));
    assert(T == length(D(subj).path));
    %T = index(3); % first half only

    init_fn = @() MH_init(D(subj), h);
    choice_fn = @(t, particle) MH_choice(t, particle, D_full(subj), h);
    update_fn = @(t, particle) MH_update(t, particle, D_full(subj), h, nsamples, T);

    clear particles;
    clear w;
    clear liks;
    clear lik;
    clear log_w;

    D(subj).tasks.s = [];
    D(subj).tasks.g = [];

    % init particles
    for i = 1:num_particles
        particles(i) = init_fn(); 
        w(i) = 1; % b/c sample from prior
    end
    w = w / sum(w);

    for t = 1:T
        % get choice probabilities
        for i = 1:num_particles
            [lik(i)] = choice_fn(t, particles(i)); % choice likelihood, P(a|D,H)
        end
        liks(t) = sum(lik .* w); % marginalize choice probability over particles

        if ismember(t, index)
            fname = sprintf('%s_subj=%d_t=%d.mat', filename, subj, t);
            fprintf('saving %s\n', fname);
            save(fname, '-v7.3');
        end

        % include observation t in data
        D(subj).tasks.s(t) = D_full(subj).tasks.s(t);
        D(subj).tasks.g(t) = D_full(subj).tasks.g(t);

        % likelihood weighing, based on observations: P(D|H)
        for i = 1:num_particles
            %log_w(i) = loglik_c(particles(i).H, D(subj), h);
            log_w(i) = pf_loglik(t, particles(i), D(subj), h);
        end
        log_w = log_w - logsumexp(log_w);
        w = exp(log_w);
        w = w / sum(w);

        particles = resample_particles(particles, w);

        % rejuvinate particles
        for i = 1:num_particles
            particles(i) = update_fn(t, particles(i));
        end
    end

    lme(subj,:) = sum(log(liks));
    lme_probes(subj,:) = sum(log(liks(index)));

    toc

    save(filename, '-v7.3');
end

save(filename, '-v7.3');
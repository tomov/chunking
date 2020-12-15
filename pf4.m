% SIR for experiment 3, like isl_MH but using observation likelihood
% c/p isl_MH.m

h = init_hyperparams;
nsamples = 1000;
num_particles = 1000;

index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..; ORDER CRUCIAL

sem = @(x) std(x) / sqrt(length(x));


% from model_all_data

filename = sprintf('mat/pf4_alpha=%.4f_nsamples=%d_div_eps=%.4f_last_np=%d.mat', h.alpha, nsamples, h.eps, num_particles);
[~, name] = system('hostname');
if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
    % local
    filename
else
    % cannon
    filename = fullfile(getenv('MY_SCRATCH'), 'chunking', filename);
    filename
end


[D, filenames] = init_Ds_from_data('exp/results/exp_v2_3_subway10_unlearn_circ', true);
D_full = D;



for subj = 1:length(D) % for each subject
    fprintf('infer H: subject %d\n', subj);

    tic

    T = length(D(subj).tasks.s);
    assert(T == 246);
    assert(T == length(D(subj).tasks.g));
    assert(T == length(D(subj).path));
    %T = index(3); % first half only

    init_fn = @() MH_init(D(subj), h);
    choice_fn = @(t, particle) pf_lik(t, particle, D(subj), h);
    update_fn = @(t, particle) MH_update(t, particle, D(subj), h, nsamples, T);

    results = forward(T, num_particles, init_fn, choice_fn, update_fn, filename, index, subj);

   % lme(subj,:) = sum(log(results.liks)); this is wrong
   blah blah
    lme_probes(subj,:) = sum(log(results.liks(index)));

    toc
    
    save(filename, '-v7.3');
end


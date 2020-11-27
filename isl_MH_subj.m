% ISL for experiment 3, for 1 subj
% c/p isl_MH.m

function isl_MH_subj(subj, num_particles)

h = init_hyperparams;
nsamples = 10000;

if ~exist('num_particles', 'var')
    num_particles = 10;
end

index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..; ORDER CRUCIAL

sem = @(x) std(x) / sqrt(length(x));


% from model_all_data

[~, name] = system('hostname');
filename = sprintf('mat/isl_MH_subj_%d_alpha=%.4f_nsamples=%d_div_eps=%.4f_last_np=%d.mat', subj, h.alpha, nsamples, h.eps, num_particles);
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


fprintf('infer H: subject %d, %s\n', subj, filenames{subj});

    tic

    T = length(D(subj).tasks.s);
    assert(T == 246);
    assert(T == length(D(subj).tasks.g));
    assert(T == length(D(subj).path));

    init_fn = @() MH_init(D(subj), h);
    choice_fn = @(t, particle) MH_choice(t, particle, D(subj), h);
    update_fn = @(t, particle) MH_update(t, particle, D(subj), h, nsamples, T);

    results(subj) = forward(T, num_particles, init_fn, choice_fn, update_fn);

    lme(subj,:) = sum(log(results(subj).liks));
    lme_probes(subj,:) = sum(log(results(subj).liks(index)));

    toc

save(filename, '-v7.3');

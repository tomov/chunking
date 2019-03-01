
clear all;
close all;
%[data, Ts] = load_data('exp/results', 246); % for exp_v2_3 (subway 10 unlearn)

% TODO cleanup first

h = init_hyperparams;
nsamples = 100;

% from model_all_data
D = init_Ds_from_data('exp/results/exp_v2_3_subway10_unlearn');
D_full = D;

% from analyze_data
% for exp_v2_3 subway 10 unlearn 
start = [6 6 6 6 6 6];
goal = [1 1 1 1 1 1];
%ordinal = [1 2 3 4 5]; <-- don't use that, e.g. if we happened to have 6 1 by chance
index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..
nexts = [
7 5;
7 5;
7 5;
7 5;
7 5;
7 5
];

for subj = 1:length(D) % for each subject
    D(subj).tasks.s = [];
    D(subj).tasks.g = [];

    [samples, post] = sample(D(subj), h, nsamples);
    H(subj) = samples(end);

    t = 1;
    for i = 1:length(index)
        while t < index(i)
            D(subj).task.s = [D(subj).task.s D_full(subj).tasks.s(t)];
            D(subj).task.g = [D(subj).task.g D_full(subj).tasks.g(t)];
            t = t + 1;
        end
        
        s = start(i);
        g = goal(i);
        assert(s == D_full(subj).tasks.s(t));
        assert(g == D_full(subj).tasks.g(t));
        
        [samples, post] = sample(D(subj), h, nsamples);
    end
end

%
% analysis reviewer suggested
%
% For example, they could make use of the variance in the random trials to perform more convincing analyses: to disprove the previous simple association of 6-5 account, they could test if there is a relationship between the number of 6-7, 7-6, 6-5, 5-6 encounters to the choice at test.

%[data, Ts] = load_data('exp/results/exp_v2_3_subway10_unlearn_circ', 246, false); % for exp_v2_3 (subway 10 unlearn)
%[data, Ts] = load_data('exp/results/subway10_repro', 83); % for subway 10 repro TODO change phase = 2!

% save analyze_reviewer.mat
load analyze_reviewer.mat

train67 = []; % # of 6-7 transitions during training for each subj
train65 = []; % # of 6-5 transitions during training for each subj

test67 = []; % whether subject went to 6-7 (as opposed to 6-5) on test trial

for subj = 1:size(data,1) % for each subject


    train67(subj,:) = 0;
    train65(subj,:) = 0;
    train76(subj,:) = 0;
    train56(subj,:) = 0;

    phase = 1; % training
    for i = 1:length(data(subj, phase).s) % for each trial 
        for j = 2:length(data(subj, phase).path{i})
            u = data(subj, phase).path{i}(j-1);
            v = data(subj, phase).path{i}(j);
            if u == 6 && v == 5
                train65(subj,:) = train65(subj,:) + 1;
            elseif u == 6 && v == 7
                train67(subj,:) = train67(subj,:) + 1;
            elseif u == 5 && v == 6
                train56(subj,:) = train56(subj,:) + 1;
            elseif u == 7 && v == 6
                train76(subj,:) = train76(subj,:) + 1;
            end
        end
    end

    phase = 2; % test
    for i = 1:length(data(subj, phase).s) % for each trial 

        if data(subj, phase).s(i) ~= 6 || data(subj, phase).g(i) ~= 1
            continue;
        end

        u = data(subj, phase).path{i}(1);
        v = data(subj, phase).path{i}(2);
        assert(u == 6);
        if v == 5
            test67(subj,:) = 0;
        elseif v == 7
            test67(subj,:) = 1;
        end
    end
end


tbl = table(test67, train67, train65, train76, train56);

formula = 'test67 ~ -1 + train76 + train56 + train67 + train65';
result = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result);

result

for i = 1:4
    H = [0 0 0 0];
    H(i) = 1;
    [p, F, DF1, DF2] = coefTest(result, H);
    fprintf('is coefficient for %d significant?  = %f, p = %f, F(%d,%d) = %f\n', i, H * beta, p, DF1, DF2, F);
end


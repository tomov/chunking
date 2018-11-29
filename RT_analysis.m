clear;


sem = @(x) std(x) / sqrt(length(x));

% --------- MAP view -------------

% for subway 10
%
action_chunk_transitions{1} = [1 2; 2 3; 4 5; 5 6; 10 9; 9 8; 8 7];
state_chunk_transitions{1} = flip(action_chunk_transitions{1},2);
bridges{1} = [3 4; 4 3; 1 10; 10 1; 7 6; 6 7];
dirname{1} = 'exp/results/subway10_map/';
nrows{1} = 81;

% for subway 9 (bad chunks)
%
action_chunk_transitions{2} = [1 2; 2 3; 4 5; 5 6; 9 8; 8 7];
state_chunk_transitions{2} = flip(action_chunk_transitions{2},2);
bridges{2} = [3 4; 4 3; 1 9; 9 1; 7 6; 6 7];
dirname{2} = 'exp/results/subway9_map/';
nrows{2} = 81;

% for subway 9 (good chunks)
%
action_chunk_transitions{3} = [1 9; 8 7; 7 6; 2 3; 3 4; 4 5];
state_chunk_transitions{3} = flip(action_chunk_transitions{3},2);
bridges{3} = [8 9; 9 8; 6 5; 5 6; 1 2; 2 1];
dirname{3} = 'exp/results/subway9_map_goodchunks/';
nrows{3} = 81;

% --------------- adjacent view ---------------

% for subway 10
%
action_chunk_transitions{4} = [1 2; 2 3; 4 5; 5 6; 10 9; 9 8; 8 7];
state_chunk_transitions{4} = flip(action_chunk_transitions{4},2);
bridges{4} = [3 4; 4 3; 1 10; 10 1; 7 6; 6 7];
dirname{4} = 'exp/results/subway10_repro';
nrows{4} = 83;

% for subway 9 (bad chunks)
%
action_chunk_transitions{5} = [1 2; 2 3; 4 5; 5 6; 9 8; 8 7];
state_chunk_transitions{5} = flip(action_chunk_transitions{5},2);
bridges{5} = [3 4; 4 3; 1 9; 9 1; 7 6; 6 7];
dirname{5} = 'exp/results/subway9/';
nrows{5} = 81;

% for subway 9 (good chunks)
%
action_chunk_transitions{6} = [1 9; 8 7; 7 6; 2 3; 3 4; 4 5];
state_chunk_transitions{6} = flip(action_chunk_transitions{6},2);
bridges{6} = [8 9; 9 8; 6 5; 5 6; 1 2; 2 1];
dirname{6} = 'exp/results/subway9_goodchunks/';
nrows{6} = 81;



action_chunk_RTs = [];
state_chunk_RTs = [];
bridge_RTs = [];

% aggregate across datasets for more power
%
for f = 1:length(dirname)
    fprintf('\n\n ---------------- Data dir %s -------------- \n\n', dirname{f});

    [data, Ts] = load_data(dirname{f}, nrows{f});

    % RT analysis
    %
    for subj = 1:size(data,1) % for each subject
        phase = 1; % training phase
        for i = 1:length(data(subj, phase).s) % for each trial
            RTs = data(subj,phase).RTs{i};
            path = data(subj,phase).path{i};
            if length(RTs) == length(path) % we log all key presses but not all of them are moves... oops
                for j = 2:length(path) - 1 % skip first RT; it's always slow
                    RT = RTs(j);
                    u = path(j);
                    v = path(j+1);
                    if any(ismember(bridges{f}, [u v], 'rows'))
                        bridge_RTs = [bridge_RTs RT];
                    elseif any(ismember(action_chunk_transitions{f}, [u v], 'rows'))
                        action_chunk_RTs = [action_chunk_RTs RT];
                    else
                        assert(any(ismember(state_chunk_transitions{f}, [u v], 'rows')));
                        state_chunk_RTs = [state_chunk_RTs RT];
                    end
                end
            else
                fprintf('skipping subj %d trial %d\n', subj, i);
            end
        end
    end
end

[h, p, ci, stats] = ttest2(action_chunk_RTs, bridge_RTs);
fprintf('bridges: %.3f +- %.3f ms\n', mean(bridge_RTs), sem(bridge_RTs));
fprintf('action chunks: %.3f +- %.3f ms\n', mean(action_chunk_RTs), sem(action_chunk_RTs));
fprintf('t(%d) = %.3f, p = %f\n', stats.df, stats.tstat, p);

fprintf('\n\n');

[h, p, ci, stats] = ttest2(state_chunk_RTs, bridge_RTs);
fprintf('bridges: %.3f +- %.3f ms\n', mean(bridge_RTs), sem(bridge_RTs));
fprintf('state chunks: %.3f +- %.3f ms\n', mean(state_chunk_RTs), sem(state_chunk_RTs));
fprintf('t(%d) = %.3f, p = %f\n', stats.df, stats.tstat, p);

figure;
m = [mean(action_chunk_RTs) mean(state_chunk_RTs) mean(bridge_RTs)]; 
se = [sem(action_chunk_RTs) sem(state_chunk_RTs) sem(bridge_RTs)]; 
bar(m);
hold on;
errorbar(m, se, 'linestyle', 'none');
ylabel('RT (ms)');
xticklabels({'action chunks', 'state chunks', 'bridges'});

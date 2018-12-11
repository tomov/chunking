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
html{1} = 'exp/exp_v3_2.html';

% for subway 9 (bad chunks)
%
action_chunk_transitions{2} = [1 2; 2 3; 4 5; 5 6; 9 8; 8 7];
state_chunk_transitions{2} = flip(action_chunk_transitions{2},2);
bridges{2} = [3 4; 4 3; 1 9; 9 1; 7 6; 6 7];
dirname{2} = 'exp/results/subway9_map/';
nrows{2} = 81;
html{2} = 'exp/exp_v3_3.html';

% for subway 9 (good chunks)
%
action_chunk_transitions{3} = [1 9; 8 7; 7 6; 2 3; 3 4; 4 5];
state_chunk_transitions{3} = flip(action_chunk_transitions{3},2);
bridges{3} = [8 9; 9 8; 6 5; 5 6; 1 2; 2 1];
dirname{3} = 'exp/results/subway9_map_goodchunks/';
nrows{3} = 81;
html{3} = 'exp/exp_v3_6.html';

% --------------- adjacent view ---------------

% for subway 10
%
action_chunk_transitions{4} = [1 2; 2 3; 4 5; 5 6; 10 9; 9 8; 8 7];
state_chunk_transitions{4} = flip(action_chunk_transitions{4},2);
bridges{4} = [3 4; 4 3; 1 10; 10 1; 7 6; 6 7];
dirname{4} = 'exp/results/subway10_repro';
nrows{4} = 83;
html{4} = 'exp/exp_v1.html';

% for subway 9 (bad chunks)
%
action_chunk_transitions{5} = [1 2; 2 3; 4 5; 5 6; 9 8; 8 7];
state_chunk_transitions{5} = flip(action_chunk_transitions{5},2);
bridges{5} = [3 4; 4 3; 1 9; 9 1; 7 6; 6 7];
dirname{5} = 'exp/results/subway9/';
nrows{5} = 81;
html{5} = 'exp/exp_v1_3.html';

% for subway 9 (good chunks)
%
action_chunk_transitions{6} = [1 9; 8 7; 7 6; 2 3; 3 4; 4 5];
state_chunk_transitions{6} = flip(action_chunk_transitions{6},2);
bridges{6} = [8 9; 9 8; 6 5; 5 6; 1 2; 2 1];
dirname{6} = 'exp/results/subway9_goodchunks/';
nrows{6} = 81;
html{6} = 'exp/exp_v1_5.html';




action_chunk_RTs = [];
state_chunk_RTs = [];
bridge_RTs = [];

% aggregate across datasets for more power
%
for f = 1:3  %length(dirname)
    fprintf('\n\n ---------------- Data dir %s -------------- \n\n', dirname{f});

    [data, Ts] = load_data(dirname{f}, nrows{f});

    % TODO this works only for the full map view; the others have weird shit like rotations
    ex_noflipped = readExp(html{f});
    ex_flipped = ex_noflipped;
    ex_flipped.adj(:,1) = ex_noflipped.adj(:,3);
    ex_flipped.adj(:,3) = ex_noflipped.adj(:,1);
    move_keys = [39 38 37 40];

    % RT analysis
    %
    for subj = 1:size(data,1) % for each subject
        phase = 1; % training phase

        flipped_count = 0;
        noflipped_count = 0;

        % first pass -- figure out which keys are movements
        % skip trial when there are more keys/RTs than button presses.... whoops
        for i = 1:length(data(subj, phase).s) % for each trial
            RTs = data(subj,phase).RTs{i};
            path = data(subj,phase).path{i};
            keys = data(subj,phase).keys{i};
            if length(RTs) == length(path) % we log all key presses but not all of them are moves... oops
                for j = 2:length(path) - 1 % skip first RT; it's always slow
                    RT = RTs(j);
                    u = path(j);
                    v = path(j+1);
                    key = keys(j);
                    dir = find(move_keys == key);

                    if ex_noflipped.adj(u, dir) ~= v % can't be non-flipped
                        flipped_count = flipped_count + 1;
%                        fprintf('impossible nonflipped: %d press %d (%d) -> %d when adj is %d\n', u, key, dir, v, ex_noflipped.adj(u, dir));
                    end
                    if ex_flipped.adj(u, dir) ~= v % can't be flipped
                        noflipped_count = noflipped_count + 1;
%                        fprintf('impossible flipped: %d press %d (%d) -> %d when adj is %d\n', u, key, dir, v, ex_flipped.adj(u, dir));
                    end

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

        % figure out if subject did flipped or non-flipped version of experiment
        assert(flipped_count == 0 || noflipped_count == 0);
        if flipped_count > 0
            ex = ex_flipped;
        else
            ex = ex_noflipped;
        end

        % second pass pass -- deal with skipped trials...
        for i = 1:length(data(subj, phase).s) % for each trial
            path = data(subj,phase).path{i};
            all_RTs = data(subj,phase).RTs{i};
            all_keys = data(subj,phase).keys{i};
            RTs = [];
            keys = [];
            if length(RTs) ~= length(path)
                k = 1;
                % remove incorrect key presses 
                for j = 1:length(path)-1
                    u = path(j);
                    while true
                        if k > length(all_keys)
                            break;
                        end
                        key = all_keys(k);
                        dir = find(move_keys == key);
                        if ~isempty(dir) && ex.adj(u, dir) >= 0
                            break;
                        end
                        k = k + 1;
                    end
                    assert(k <= length(all_keys));
                    keys = [keys all_keys(k)];
                    RTs = [RTs all_RTs(k)];
                end
                % special care for last key press (space)
                keys = [keys all_keys(end)];
                RTs = [RTs all_RTs(end)];
                assert(length(RTs) == length(path));
                assert(length(keys) == length(path));

                % NOW do the RT stuff like in the first pass...
                for j = 2:length(path) - 1 % skip first RT; it's always slow
                    RT = RTs(j);
                    u = path(j);
                    v = path(j+1);
                    key = keys(j);

                    if any(ismember(bridges{f}, [u v], 'rows'))
                        bridge_RTs = [bridge_RTs RT];
                    elseif any(ismember(action_chunk_transitions{f}, [u v], 'rows'))
                        action_chunk_RTs = [action_chunk_RTs RT];
                    else
                        assert(any(ismember(state_chunk_transitions{f}, [u v], 'rows')));
                        state_chunk_RTs = [state_chunk_RTs RT];
                    end
                end
            end
        end
    end
end

%{
[h, p, ci, stats] = ttest2(action_chunk_RTs, bridge_RTs);
fprintf('bridges: %.3f +- %.3f ms\n', mean(bridge_RTs), sem(bridge_RTs));
fprintf('vs. action chunks: %.3f +- %.3f ms\n', mean(action_chunk_RTs), sem(action_chunk_RTs));
fprintf('t(%d) = %.3f, p = %f\n', stats.df, stats.tstat, p);

fprintf('\n\n');

[h, p, ci, stats] = ttest2(state_chunk_RTs, bridge_RTs);
fprintf('bridges: %.3f +- %.3f ms\n', mean(bridge_RTs), sem(bridge_RTs));
fprintf('vs. state chunks: %.3f +- %.3f ms\n', mean(state_chunk_RTs), sem(state_chunk_RTs));
fprintf('t(%d) = %.3f, p = %f\n', stats.df, stats.tstat, p);

fprintf('\n\n');

[h, p, ci, stats] = ttest2(state_chunk_RTs, action_chunk_RTs);
fprintf('action_chunks: %.3f +- %.3f ms\n', mean(action_chunk_RTs), sem(action_chunk_RTs));
fprintf('vs. state chunks: %.3f +- %.3f ms\n', mean(state_chunk_RTs), sem(state_chunk_RTs));
fprintf('t(%d) = %.3f, p = %f\n', stats.df, stats.tstat, p);
%}

RTs = [action_chunk_RTs state_chunk_RTs bridge_RTs];
[group{1:length(action_chunk_RTs)}] = deal('action chunks');
[group{length(action_chunk_RTs)+1:length(action_chunk_RTs)+length(state_chunk_RTs)}] = deal('state chunks');
[group{length(action_chunk_RTs)+length(state_chunk_RTs)+1:length(RTs)}] = deal('bridges');
[p, tbl, stats] = anova1(RTs, group);

[c,~,~,names] = multcompare(stats);
res = [{'group 1', 'group 2', 'lower CI', 'mean diff', 'upper CI', 'p-value'}; names(c(:,1)), names(c(:,2)), num2cell(c(:,3:6))];
res


figure;
m = [mean(action_chunk_RTs) mean(state_chunk_RTs) mean(bridge_RTs)]; 
se = [sem(action_chunk_RTs) sem(state_chunk_RTs) sem(bridge_RTs)]; 
bar(m);
hold on;
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
ylabel('RT (ms)');
xticklabels({'action chunks', 'state chunks', 'bridges'});

save('RT_analysis.mat');


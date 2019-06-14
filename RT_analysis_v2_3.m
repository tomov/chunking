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


% for subway_unlearn 
%
action_chunk_transitions{7} = [2 3; 4 5; 10 9; 9 8; 8 7];
state_chunk_transitions{7} = flip(action_chunk_transitions{1},2);
bridges{7} = [1 2; 2 1; 3 4; 4 3; 5 6; 6 5; 1 10; 10 1; 7 6; 6 7];
dirname{7} = 'exp/results/exp_v2_3_subway10_unlearn_circ';
nrows{7} = 246;
html{7} = 'exp/exp_v2_3.html';




action_chunk_RTs = [];
state_chunk_RTs = [];
bridge_RTs = [];

% for fitglme
rt = []; % all RTs
type = []; % RT category 1 = action chunk, 2 = state chunk, 3 = bridge
subject = [];  % subject
experiment = []; %

for f = 7:7  %length(dirname)
    fprintf('\n\n ---------------- Data dir %s -------------- \n\n', dirname{f});

    %[data, Ts] = load_data(dirname{f}, nrows{f});
    %save RT_analysis_tmp.mat
    load RT_analysis_tmp.mat


    move_keys = [39 38 37 40];

    % RT analysis
    %
    for subj = 1:size(data,1) % for each subject
        phase = 1; % training phase

        ex_ruled_out = [0 0 0 0]; % rule out different orientations in exs{} one by one, i.e. see if any moves rule those out

        % first pass -- figure out which keys are movements
        % skip trial when there are more keys/RTs than button presses.... whoops
        for i = 1:length(data(subj, phase).s) % for each trial
            RTs = data(subj,phase).RTs{i};
            path = data(subj,phase).path{i};
            keys = data(subj,phase).keys{i};
            valid_keys = data(subj,phase).valid_keys{i};

            % only use valid keys that resulted in movements
            RTs = [RTs(valid_keys) RTs(end)];
            keys = [keys(valid_keys) keys(end)]; % last key is always space

            assert(length(RTs) == length(path));
            for j = 2:length(path) - 1 % skip first RT; it's always slow; also skip last one, it's always space
                RT = RTs(j);
                u = path(j);
                v = path(j+1);
                key = keys(j);
                dir = find(move_keys == key);

                % for fitglme
                rt = [rt; RT];
                subject = [subject; subj];
                experiment = [experiment; f];

                if any(ismember(bridges{f}, [u v], 'rows'))
                    bridge_RTs = [bridge_RTs RT];
                    type = [type; 3];
                elseif any(ismember(action_chunk_transitions{f}, [u v], 'rows'))
                    action_chunk_RTs = [action_chunk_RTs RT];
                    type = [type; 1];
                else
                    assert(any(ismember(state_chunk_transitions{f}, [u v], 'rows')));
                    state_chunk_RTs = [state_chunk_RTs RT];
                    type = [type; 2];
                end
            end
        end

    end
end

type = categorical(type); % important! for fitglme
experiment = categorical(experiment);

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

save('RT_analysis_forglme_v2_3.mat');
% see fig_RT_analysis_2.m


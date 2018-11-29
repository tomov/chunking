% call analyze_data.m first

% for subway 10
%
action_chunk_transitions = [
1 2;
2 3;
4 5;
5 6;
10 9;
9 8;
8 7];

bridges = [
3 4;
4 3;
1 10;
10 1;
7 6;
6 7];

state_chunk_transitions = flip(action_chunk_transitions,2);

% RT analysis
%
action_chunk_RTs = [];
state_chunk_RTs = [];
bridge_RTs = [];
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
                if any(ismember(bridges, [u v], 'rows'))
                    bridge_RTs = [bridge_RTs RT];
                elseif any(ismember(action_chunk_transitions, [u v], 'rows'))
                    action_chunk_RTs = [action_chunk_RTs RT];
                else
                    assert(any(ismember(state_chunk_transitions, [u v], 'rows')));
                    state_chunk_RTs = [state_chunk_RTs RT];
                end
            end
        else
            fprintf('skipping subj %d trial %d\n', subj, i);
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

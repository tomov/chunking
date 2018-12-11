
load('RT_analysis.mat');

RTs = [action_chunk_RTs state_chunk_RTs bridge_RTs];
[group{1:length(action_chunk_RTs)}] = deal('action chunks');
[group{length(action_chunk_RTs)+1:length(action_chunk_RTs)+length(state_chunk_RTs)}] = deal('state chunks');
[group{length(action_chunk_RTs)+length(state_chunk_RTs)+1:length(RTs)}] = deal('bridges');
[p, tbl, stats] = anova1(RTs, group);



figure;
m = [mean(action_chunk_RTs) mean(state_chunk_RTs) mean(bridge_RTs)]; 
se = [sem(action_chunk_RTs) sem(state_chunk_RTs) sem(bridge_RTs)]; 
bar(m);
hold on;
errorbar(m, se, 'linestyle', 'none', 'color', 'black');
ylabel('RT (ms)');
xticklabels({'action chunk', 'state chunk', 'boundary'});


print('figures/RTs.pdf', '-dpdf');



[c,~,~,names] = multcompare(stats);
res = [{'group 1', 'group 2', 'lower CI', 'mean diff', 'upper CI', 'p-value'}; names(c(:,1)), names(c(:,2)), num2cell(c(:,3:6))];
res


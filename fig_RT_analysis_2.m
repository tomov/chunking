% using fitglme mixed effects

load('RT_analysis_forglme.mat');

type = categorical(type); % 1 = action chunk, 2 = state chunk, 3 = bridge
experiment = categorical(experiment);

tbl = table(rt, type, subject, experiment);

which = rt > 0;
tbl = tbl(which,:);
tbl.rt = log(tbl.rt); % log transform

formula = 'rt ~ 1 + type + (1 + type | subject)';

result = fitglme(tbl, formula, 'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace');
[beta, names, stats] = fixedEffects(result);

H = [0 -1 1 ]; % bridge - state chunk
[p, F, DF1, DF2] = coefTest(result, H);
fprintf('fitglme bridge - state chunk contrast: = %f (expect positive), p = %f, F(%d,%d) = %f\n', H * beta, p, DF1, DF2, F);




figure;

b = rt(type == categorical(3) & rt > 0);
s = rt(type == categorical(2) & rt > 0);

histogram(log(s));
histogram(log(b));
legend({'state chunk', 'bridge'});
ylabel('# presses');
xlabel('log RT (log ms)');
%figure;
%m = [mean(action_chunk_RTs) mean(state_chunk_RTs) mean(bridge_RTs)]; 
%se = [sem(action_chunk_RTs) sem(state_chunk_RTs) sem(bridge_RTs)]; 
%bar(m);
%hold on;
%errorbar(m, se, 'linestyle', 'none', 'color', 'black');
%ylabel('RT (ms)');
%xticklabels({'action chunk', 'state chunk', 'boundary'});



clear all;
close all;


prefix = 'pf2'; % isl_MH or pf2
%prefix = 'isl_MH'; % isl_MH or pf2

template = sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=103.mat', prefix);
filename = sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', prefix);

[~, name] = system('hostname');
if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
    % local
    template
    filename

else
    % cannon
    template = fullfile(getenv('MY_SCRATCH'), 'chunking', template);
    filename = fullfile(getenv('MY_SCRATCH'), 'chunking', filename);

    template
    filename

end

files = dir(template);

load(filename, 'D');


H = [1 1 1 2 2 2 3 3 3 3]; % ground truth

for i = 1:length(files)
    files(i).name

    load(fullfile(files(i).folder, files(i).name), 'subj', 'particles');

    for j = 1:length(particles)
        rs(j) = rand_index(H, particles(j).H.c, 'adjusted');
    end
    r(subj) = mean(rs); % mean rand index
    c(subj) = D(subj).path{103}(2); % choice on last probe trial in first half
end

%save('isl_fig_isl_MH_adj.mat')
save('isl_fig_pf2_adj.mat')

%load('mat/isl_fig_isl_MH.mat');
%load('mat/isl_fig_pf2.mat');

sem = @(x) std(x) / sqrt(length(x));

ms = [mean(r(c == 5)) mean(r(c == 7))];
sems = [sem(r(c == 5)) sem(r(c == 7))];

figure;

hold on;
bar([1 2], ms);
errorbar([1 2], ms, sems, 'linestyle', 'none', 'color', 'black');

x(c == 5) = 1;
x(c == 7) = 2;
swarmchart(x, r);

xlabel('subject choice on trial 103');
ylabel('avg particle rand index');
set(gca, 'xtick', [1 2]);
xticklabels({'5', '7'});


%[h,p,ci,stats] = ttest2(r(c == 5), r(c == 7))

[p,h,stats] = ranksum(r(c == 5), r(c == 7))


clear all;
close all;


prefix = 'pf2'; % isl_MH or pf2
cachefile = 'mat/isl_fig_pf2_6.mat';

%{
H = [1 1 1 2 2 2 3 3 3 3]; % ground truth

index = [34 68 103 47+103 94+103 143+103]; % from html -- @ ..; ORDER CRUCIAL

templates = {sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=34.mat', prefix), 
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=68.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=103.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=150.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=197.mat', prefix),
             sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat_subj=*_t=246.mat', prefix)};
filename = sprintf('mat/%s_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', prefix);

[~, name] = system('hostname');
if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
    % local
    templates
    filename

else
    % cannon
    for t = 1:length(templates)
        templates{t} = fullfile(getenv('MY_SCRATCH'), 'chunking', templates{t});
    end
    filename = fullfile(getenv('MY_SCRATCH'), 'chunking', filename);

    templates
    filename

end

for t = 1:length(templates)

    files = dir(templates{t});

    load(filename, 'D');

    for i = 1:length(files)
        files(i).name

        load(fullfile(files(i).folder, files(i).name), 'subj', 'particles');

        rr = nan(size(particles));
        for j = 1:length(particles)
            rr(j) = rand_index(H, particles(j).H.c, 'adjusted');
        end
        rs{t}(subj) = mean(rr); % mean rand index
        cs{t}(subj) = D(subj).path{index(t)}(2); % choice on last probe trial in first half
    end
end

cachefile
save(cachefile);

%}
load(cachefile);

%r = [rs{1} rs{2} rs{3}];
%c = [cs{1} cs{2} cs{3}];

sem = @(x) std(x) / sqrt(length(x));

figure;

% ------------ 1st half

subplot(1,2,1);

r = rs{3};
c = cs{3};

ms = [mean(r(c == 5)) mean(r(c == 7))];
sems = [sem(r(c == 5)) sem(r(c == 7))];

bar([1 2], ms);
hold on;
errorbar([1 2], ms, sems, 'linestyle', 'none', 'color', 'black');
x(c == 5) = 1;
x(c == 7) = 2;
swarmchart(x, r);

title('first half');
xlabel('subject choice on trial 103');
ylabel('avg particle rand index');
set(gca, 'xtick', [1 2]);
xticklabels({'5', '7'});

[h,p,ci,stats] = ttest2(r(c == 5), r(c == 7))
[p,h,stats] = ranksum(r(c == 5), r(c == 7))


% ------------ 2nd half


subplot(1,2,2);

r = rs{6};
c = cs{6};

ms = [mean(r(c == 5)) mean(r(c == 7))];
sems = [sem(r(c == 5)) sem(r(c == 7))];

bar([1 2], ms);
hold on;
errorbar([1 2], ms, sems, 'linestyle', 'none', 'color', 'black');
x(c == 5) = 1;
x(c == 7) = 2;
swarmchart(x, r);

title('second half');
xlabel('subject choice on trial 246');
ylabel('avg particle rand index');
set(gca, 'xtick', [1 2]);
xticklabels({'5', '7'});

[h,p,ci,stats] = ttest2(r(c == 5), r(c == 7))
[p,h,stats] = ranksum(r(c == 5), r(c == 7))


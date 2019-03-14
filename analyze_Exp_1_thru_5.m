clear all;
close all;

init_all_plots;

filename = 'analyze_Exp_1_thru_5.mat';
filename

sem = @(x) std(x) / sqrt(length(x));

for i = 1:length(pl)
    for j = 1:length(pl(i).dirnames)
        if ~isnan(pl(i).m(j))
            continue;
        end

        [data, Ts] = load_data(pl(i).dirnames{j}, pl(i).nrows(j));
        %pl(i).data{j} = data;
        %pl(i).Ts{j} = Ts;

        % TODO dedupe w/ analyze_data.m
        s = [];
        g = [];
        len = [];
        group = [];
        dir = []; % direction = 2nd state on path
        ord = []; % ordinal of trial type within phase (e.g. "first 1->6", "second 1->6", etc)
        subj_group = [];
        subj_len = [];
        s_id = [];
        for subj = 1:size(data,1)
            phase = 2;
            for k = 1:length(data(subj, phase).s)
                which = find(data(subj, phase).s == data(subj, phase).s(k) & data(subj, phase).g == data(subj, phase).g(k));
                clear o;
                o(which) = find(which);
                ord = [ord; o(k)];
                s = [s; data(subj, phase).s(k)];
                g = [g; data(subj, phase).g(k)];
                len = [len; data(subj, phase).len(k)];
                dir = [dir; data(subj, phase).path{k}(2)];
                group = [group; data(subj, phase).group(k)];
                s_id = [s_id; subj];
            end
            subj_group = [subj_group; data(subj,1).group(1)];
            subj_len = [subj_len; mean(data(subj, 1).len)];
        end
      
        % TODO dedupe w/ analyze_data.m
        which = s == pl(i).starts(j);
        move = dir(which);
        c1 = sum(move == pl(i).nexts(j)); % count 1
        c2 = sum(move ~= pl(i).nexts(j)); % count 2
        n = sum(which);
        switch pl(i).tests(j)
            case 1 % right-tailed
                p = 1 - binocdf(c1, n, 0.5);
                assert(false);
            case 2 % left-tailed
                p = binocdf(c1, n, 0.5);
                assert(false);
            case 3 % two-tailed
                p = 2 * binocdf(min(c1,c2), n, 0.5);
            otherwise
                assert(false);
        end

        y = binoinv([0.025 0.975], n, 0.5);
        pl(i).ci(j) = (y(2) - y(1)) / 2;
        pl(i).n(j) = n;
        pl(i).m(j) = c1;
        pl(i).p(j) = p;
    end
end

clear Ts;
clear data;
save(filename, '-v7.3');

%load('analyze_all_data.mat');

plot_all_data;

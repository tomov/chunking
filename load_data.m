function [data, Ts] = load_data

    dirname = 'exp/results'; 
    bad_dirname = 'exp/results/bad';

    expected_number_of_rows = 83;
    %dirname = 'exp/results/subway10_repro';  % expected = 83
    %dirname = 'exp/results/subway10'; % expected # rows = 100
    %dirname = 'exp/results/subway10_noadj'; % expected # rows = 110
    %dirname = 'exp/results/subway_10_randsg_WRONG'; % expected # rows = 116

    %dirname = 'exp/results/ARCHIVE/subway_10_noadj_batch2'; % expected # rows = 110
    %dirname = 'exp/results/ARCHIVE/subway10_batch2'; % expected # rows = 110

    files = dir(dirname);
    subj = 1;
    for idx = 1:length(files)
        if ~endsWith(files(idx).name, 'csv')
            continue;
        end

        filepath = fullfile(dirname, files(idx).name);
        try
            T = readtable(filepath);
        catch
            fprintf('Error reading file %s\n', files(idx).name);
            if exist('bad_dirname', 'var')
                movefile(filepath, bad_dirname);
            end
            continue;
        end
        if size(T, 1) ~= expected_number_of_rows
            fprintf('Skipping %s: it has only %d rows\n', files(idx).name, size(T,1));
            if exist('bad_dirname', 'var')
                movefile(filepath, bad_dirname);
            end
            continue;
        end
        Ts{subj} = T;

        skip_subj = false;

        % TODO dedupe with init_D_from_csv.m
        max_RT = 0;
        phase = 1;
        j = 1; % idx within phase
        for i = 1:size(T,1)
            stage = strip(T.stage{i});
            switch phase
                case 1
                    if strcmp(stage, 'test')
                        phase = 2;
                        j = 1;
                    end
                case 2
                    if strcmp(stage, 'training')
                        phase = 3;
                        j = 1;
                    end
                case 3
                    if strcmp(stage, 'test')
                        phase = 4;
                        j = 1;
                    end
            end

            s = T.start(i);
            if iscell(s)
                s = str2num(s{1});
            end
            if isempty(s)
                s = 0;
            end
            g = T.goal(i);
            if iscell(g)
                g = str2num(g{1});
            end
            if isempty(g)
                g = 0;
            end
            RTs = str2num(T.RTs{i});
            max_RT = max(max_RT, max(RTs));
            path = str2num(T.path{i});
            assert(length(path) == T.length(i));
            group = strip(T.group{i});
            RT_tot = T.RT_tot(i);
            switch group
                case 'A'
                    group = 1;
                case 'B'
                    group = 2;
                otherwise
                    assert(false);
            end
            id = T.subj_id(i);

            % skip subjects with unrealistically long paths
            %{
            if length(path) > 25
                fprintf('Skipping %s: trial %d has path length %d\n', files(idx).name, i, length(path));
                skip_subj = true;
                break;
            end
            %}

            data(subj, phase).s(j) = s;
            data(subj, phase).g(j) = g;
            data(subj, phase).path{j} = path;
            data(subj, phase).len(j) = length(path);
            data(subj, phase).group(j) = group;
            data(subj, phase).id = id;
            data(subj, phase).RTs{j} = RTs;
            data(subj, phase).RT_tot(j) = RT_tot;

            j = j + 1;
        end

        %{
        % skip subjects that didn't improve over time
        if ~skip_subj
            l = data(subj,1).len;
            first = l(1:round(length(l) * 0.10));
            last = l(end-round(length(l) * 0.10):end);
            [h, p, ci, stat] = ttest2(first, last, 'tail', 'right');
            if p > 0.1
                fprintf('Skipping %s: no improvement in path length (p = %.3f, first 20 = %.2f, last 20 = %.2f)\n', files(idx).name, p, mean(first), mean(last));
                skip_subj = true;
            end
        end
        %}

        fprintf('         max RT = %.2f s, total RT = %.2f min\n', max_RT / 1000, sum(T.RT_tot) / 1000 / 60);

        if ismember('timestamp', T.Properties.VariableNames)
            dur = T.timestamp(end) - T.timestamp(1);
            fprintf('             duration = %.2f mins\n', dur / 60);
        end

        if ~skip_subj
            subj = subj + 1;
        else
            if size(data,1) >= subj
                data(subj,:) = [];
            end
        end
    end

    save('data.mat', 'data', 'Ts');

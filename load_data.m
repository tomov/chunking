function [data, Ts] = load_data

    dirname = 'exp/results';

    files = dir(dirname);
    subj = 1;
    for i = 1:length(files)
        if ~endsWith(files(i).name, 'csv')
            continue;
        end

        try
            T = readtable(fullfile(dirname, files(i).name));
        catch
            fprintf('Error reading file %s\n', files(i).name);
            continue;
        end
        if size(T, 1) ~= 100
            fprintf('Skipping %s: it has only %d rows\n', files(i).name, size(T,1));
            continue;
        end
        Ts{subj} = T;

        % TODO dedupe with init_D_from_csv.m
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
            path = str2num(T.path{i});
            assert(length(path) == T.length(i));
            group = strip(T.group{i});
            switch group
                case 'A'
                    group = 1;
                case 'B'
                    group = 2;
                otherwise
                    assert(false);
            end
            id = T.subj_id(i);

            data(subj, phase).s(j) = s;
            data(subj, phase).g(j) = g;
            data(subj, phase).path{j} = path;
            data(subj, phase).len(j) = length(path);
            data(subj, phase).group(j) = group;
            data(subj, phase).id = id;

            j = j + 1;
        end

        subj = subj + 1;
    end

    save('data.mat', 'data', 'Ts');

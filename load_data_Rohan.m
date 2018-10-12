function [data, Ts] = load_data

    dirname = 'occluder3D_data';

    files = dir(dirname);
    subj = 1;
    for i = 1:length(files)
        if ~endsWith(files(i).name, 'csv')
            continue;
        end

        T = readtable(fullfile(dirname, files(i).name));
        Ts{subj} = T;

        % TODO dedupe with init_D_from_csv.m
        phase = 1;
        j = 1; % idx within phase
        for i = 1:size(T,1)
            type = strip(T.Type{i});
            switch phase
                case 1
                    if strcmp(type, 'test')
                        phase = 2;
                        j = 1;
                    end
                case 2
                    if strcmp(type, 'training')
                        phase = 3;
                        j = 1;
                    end
                case 3
                    if strcmp(type, 'test')
                        phase = 4;
                        j = 1;
                    end
            end

            s = T.Displayed(i);
            if iscell(s)
                s = str2num(s{1});
            end
            if isempty(s)
                s = 0;
            end
            g = T.Goals(i);
            if iscell(g)
                g = str2num(g{1});
            end
            if isempty(g)
                g = 0;
            end
            path = str2num(T.Path{i});
            group = strip(T.Group{i});
            switch group
                case 'Group A'
                    group = 1;
                case 'Group B'
                    group = 2;
                otherwise
                    assert(false);
            end
            id = T.ID(i);

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

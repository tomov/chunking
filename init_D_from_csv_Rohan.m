function D = init_D_from_csv(filename)
    assert(endsWith(filename, '.csv'));

    D.G.N = 0;
    D.tasks.s = [];
    D.tasks.g = [];

    T = readtable(filename);

    D.name = [strip(T.Group{1}), ' ', num2str(T.ID(1))];

    phase = 'training_1';
    for i = 1:size(T,1)
        type = strip(T.Type{i});
        switch phase
            case 'training_1'
                if strcmp(type, 'test')
                    phase = 'test_1';
                end
            case 'test_1'
                if strcmp(type, 'training')
                    phase = 'training_2';
                end
            case 'training_2'
                if strcmp(type, 'test')
                    phase = 'test_2';
                end
        end

        if ~strcmp(phase, 'training_1')
            % TODO 2nd half
            break;
        end

        s = T.Displayed(i);
        g = T.Goals(i);
        if iscell(s)
            s = str2num(s{1});
        end
        if iscell(g)
            g = str2num(g{1});
        end
        D.tasks.s = [D.tasks.s s];
        D.tasks.g = [D.tasks.g g];
        D.G.N = max([D.G.N s g]);

        path = strsplit(strip(T.Path{i}), ' ');
        for k = 1:length(path) - 1
            i = str2num(path{k});
            j = str2num(path{k + 1});
            D.G.E(i,j) = 1;
            D.G.E(j,i) = 1;
            D.G.N = max([D.G.N i j]);
        end
    end
end

function D = init_D_from_csv(filename, test_phase_too)
    assert(endsWith(filename, '.csv'));

    if ~exist('test_phase_too', 'var')
        test_phase_too = false;
    end

    D.G.N = 0;
    D.tasks.s = [];
    D.tasks.g = [];

    T = readtable(filename, 'Delimiter', ',');

    D.name = [strip(T.group{1}), ' ', num2str(T.subj_id(1))];
    %D.name = num2str(T.subj_id(1));

    phase = 'training_1';
    for i = 1:size(T,1)
        stage = strip(T.stage{i});
        switch phase
            case 'training_1'
                if strcmp(stage, 'test')
                    phase = 'test_1';
                end
            case 'test_1'
                if strcmp(stage, 'train')
                    phase = 'training_2';
                end
            case 'training_2'
                if strcmp(stage, 'test')
                    phase = 'test_2';
                end
        end

        if ~strcmp(phase, 'training_1') 
            if ~strcmp(phase, 'test_1') && ~test_phase_too
                % TODO this is a hack for model_exp_v2_3
                break;
            end
            % TODO 2nd half
        end

        s = T.start(i);
        g = T.goal(i);
        if iscell(s)
            s = str2num(s{1});
        end
        if iscell(g)
            g = str2num(g{1});
        end
        D.tasks.s = [D.tasks.s s];
        D.tasks.g = [D.tasks.g g];
        D.G.N = max([D.G.N s g]);

        path = strsplit(strip(T.path{i}), ' ');
        for k = 1:length(path) - 1
            i = str2num(path{k});
            j = str2num(path{k + 1});
            D.G.E(i,j) = 1;
            D.G.E(j,i) = 1;
            D.G.N = max([D.G.N i j]);
        end
    end
end

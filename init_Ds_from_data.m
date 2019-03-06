function [D, filenames] = init_Ds_from_data(dirname, test_phase_too)


    if ~exist('test_phase_too', 'var')
        test_phase_too = false;
    end

    files = dir(dirname);
    idx = 1;
    for i = 1:length(files)
        if ~endsWith(files(i).name, 'csv')
            continue;
        end

        D(idx) = init_D_from_csv(fullfile(dirname, files(i).name), test_phase_too);
        filenames{idx} = fullfile(dirname, files(i).name);
        idx = idx + 1;
    end

end

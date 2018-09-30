function D = init_Ds_from_data(dirname)

    files = dir(dirname);
    idx = 1;
    for i = 1:length(files)
        if ~endsWith(files(i).name, 'csv')
            continue;
        end

        D(idx) = init_D_from_csv(fullfile(dirname, files(i).name));
        idx = idx + 1;
    end

end

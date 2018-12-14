function D = init_D_from_txt(filename)
    f = fopen(filename, 'r');

    D.name = fgets(f);
    D.name = strip(D.name);
    A = freadline(f, '%d %d');
    N = A(1); M = A(2);
    D.G.N = N;
    D.G.E = zeros(N, N); % TODO sparse?
    for k = 1:M
        A = freadline(f, '%d %d');
        i = A(1); j = A(2);
        D.G.E(i,j) = 1;
        D.G.E(j,i) = 1;
        D.G.edges(k,:) = [i,j];
    end

    % tasks
    %
    D.tasks.s = [];
    D.tasks.g = [];
    n = freadline(f, '%d');
    for t = 1:n
        A = freadline(f, '%d %d');
        s = A(1); g = A(2);

        if s <= 0
            s = randi(D.G.N);
        end
        while g <= 0 || g == s
            g = randi(D.G.N);
        end

        D.tasks.s = [D.tasks.s s];
        D.tasks.g = [D.tasks.g g];
    end

    % rewards
    %
    A = freadline(f, '%d');
    O = A(1);
    for i = 1:D.G.N
        D.r{i} = [];
    end
    
    for o = 1:O
        A = freadline(f, '%d %d');
        i = A(1);
        D.r{i} = [D.r{i} A(2)];
    end

    fclose(f);
end

function A = freadline(f, spec)
    line = fgets(f);
    A = sscanf(line, spec);
end


function D = init_D_from_txt(filename)
    f = fopen(filename, 'r');

    D.name = fgets(f);
    D.name = strip(D.name);
    A = freadline(f, '%d %d');
    N = A(1); M = A(2);
    D.G.N = N;
    D.G.E = zeros(N, N); % TODO sparse?
    D.G.edges = [];
    D.G.hidden_E = zeros(N, N); % TODO sparse?
    D.G.hidden_edges = [];

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

    % hidden edges
    %
    A = freadline(f, '%d');
    if ~isempty(A)
        U = A(1);
        for u = 1:U
            A = freadline(f, '%d %d');
            i = A(1);
            j = A(2);
            D.G.hidden_edges = [D.G.hidden_edges; i j];
            D.G.hidden_E(i,j) = 1;
            D.G.hidden_E(j,i) = 1;
        end
    end

    fclose(f);
end

function A = freadline(f, spec)
    line = fgets(f);
    if line == -1
        A = [];
    else
        A = sscanf(line, spec);
    end
end


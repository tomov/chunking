function D = init_D_from_txt(filename)
    f = fopen(filename, 'r');

    D.name = fgets(f);
    D.name = strip(D.name);
    A = freadline(f, '%d %d');
    N = A(1); M = A(2);
    D.G.N = N;
    D.G.E = zeros(N, N); % TODO sparse?
    D.G.I = zeros(N, N);
    for k = 1:M
        A = freadline(f, '%d %d %d');
        i = A(1); j = A(2); exists = A(3);
        D.G.E(i,j) = exists;
        D.G.E(j,i) = exists;
        D.G.I(i,j) = 1;
        D.G.I(j,i) = 1;
    end

    D.tasks.s = [];
    D.tasks.g = [];
    n = freadline(f, '%d');
    for t = 1:n
        A = freadline(f, '%d %d');
        s = A(1); g = A(2);
        D.tasks.s = [D.tasks.s s];
        D.tasks.g = [D.tasks.g g];
    end
    
    % updates
    A = freadline(f, '%d');
    num_updates = A(1);
    D.updates = [];
    for k = 1:num_updates
        A = freadline(f, '%d %d %d');
        i = A(1); j = A(2); exists = A(3);
        D.updates = [D.updates; i j exists];
    end

    fclose(f);
end

function A = freadline(f, spec)
    line = fgets(f);
    A = sscanf(line, spec);
end


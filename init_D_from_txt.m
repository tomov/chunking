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
    end
    
    for i = 1:D.G.N
        A = freadline(f, '%d');
        r = A(1);
        D.r(i) = r;
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

    fclose(f);
end

function A = freadline(f, spec)
    line = fgets(f);
    A = sscanf(line, spec);
end


function D = init_D()
    % TODO read from file

    D.G.N = 6;
    D.G.V = 1:D.G.N;
    D.G.E = zeros(D.G.N, D.G.N);

    edges = [1 3; 1 2; 2 3; 3 4; 4 5; 4 6; 5 6];
    for i = 1:size(edges, 1)
        D.G.E(edges(i,1), edges(i,2)) = 1;
        D.G.E(edges(i,2), edges(i,1)) = 1;
    end

    % TODO tasks
    D.tasks.s = [];
    D.tasks.g = [];
end

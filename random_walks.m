function paths = random_walks(D, nwalks, steps)

    paths = {};
    for i = 1:nwalks
        s = randi(D.G.N);
        path = random_walk(s, D, steps);
        paths = [paths; {path}];
    end
end

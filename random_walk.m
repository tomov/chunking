function path = random_walk(s,D,steps)

    path = [s];
    for i = 1:steps
        next = find(D.G.E(s,:));
        s = datasample(next,1);
        path = [path s];
    end
end

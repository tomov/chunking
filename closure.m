function A = closure(E)
    % compute the transitive closure of the graph G
    %
    A = (eye(size(E)) + E) ^ size(E,1) > 0;
end

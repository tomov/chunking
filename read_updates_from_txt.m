function [num_updates, V_updates, E_updates] = read_updates_from_txt(filename)
    f = fopen(filename, 'r');
    
    A = freadline(f, '%d');
    num_updates = A(1);
    V_updates = zeros(num_updates, 1);
    E_updates = {};
    
    for i = 1:num_updates
        A = freadline(f, '%d');
        num_vertices = A(1);
        V_updates(i) = num_vertices;
        A = freadline(f, '%d');
        num_edges = A(1);
        edges = {};
        for j = 1:num_edges
            A = freadline(f, '%d %d');
            u = A(1); v = A(2);
            edges{j} = [u v];
        end
        E_updates{i} = edges;
    end

    fclose(f);
end

function A = freadline(f, spec)
    line = fgets(f);
    A = sscanf(line, spec);
end

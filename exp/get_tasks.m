function get_tasks(s, g, n)
    for i = 1:length(s)
        fprintf('%d ', s(i));
    end
    fprintf('-> ');
    for i = 1:length(g)
        fprintf('%d ', g(i));
    end
    fprintf('x %d\n', n);
end

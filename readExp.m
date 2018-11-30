function ex = readExp(filename)

% see readExp() in exp_v1.js
%

fid = fopen(filename);

state = 0;
l = 1;
lines = {};
while true
    line = fgetl(fid);
    if line == -1
        break
    end

    if ~isempty(strfind(line, '</textarea>')) && state == 1
        state = 2;
    end

    if state == 1
        lines = [lines; {line}];
        disp(line);
    end

    if ~isempty(strfind(line, '<textarea id="experiment"')) && state == 0
        state = 1;
    end
end


ex.N = sscanf(lines{1}, '%d');

for i = 1:ex.N
    A = sscanf(lines{i + 1}, '%d %d %d %d');
    ex.adj(i,:) = A(1:4);
end

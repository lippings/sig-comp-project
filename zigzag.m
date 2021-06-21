function [arr] = zigzag(mat)
% Unroll 2-D matrix in zig-zag order
N = numel(mat);
[h, w] = size(mat);
arr = zeros(1, N);

level = 1;
ind = 1;
for i = 1:(w+h-1)
    if level == 1
        arr(ind) = mat(1, 1);
        ind = ind + 1;
    else
        min_row = max(1, level-w+1);
        max_row = min(level, h);

        min_col = max(1, level-h+1);
        max_col = min(level, w);
        if mod(level, 2) == 0
            row_iter = max_row:-1:min_row;
            col_iter = min_col:1:max_col;
        else
            row_iter = min_row:1:max_row;
            col_iter = max_col:-1:min_col;
        end
        
        if length(row_iter) ~= length(col_iter)
            level
        end

        for j = 1:length(row_iter)
            row = row_iter(j);
            col = col_iter(j);

            arr(ind) = mat(row, col);
            ind = ind + 1;
        end
    end
    level = level + 1;
end
end


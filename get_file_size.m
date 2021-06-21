function [size_in_bytes] = get_file_size(path)
% Get file size in bytes
s = dir(path);
size_in_bytes = s.bytes;
end


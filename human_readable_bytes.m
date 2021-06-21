function [hr_size, unit] = human_readable_bytes(bytes)
% Convert bytes into human readable form
prefixes = ["k", "M", "G", "T", "P", "Y", "Z"];
prefix = '';

for i = 1:length(prefixes)
    if bytes > 1024
        bytes = bytes / 1024;
        prefix = char(prefixes(i));
    else
        break;
    end
end

hr_size = bytes;
unit = [prefix 'B'];
end
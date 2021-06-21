function [hr_size] = human_readable_file_size(filepath, rounding)
if nargin < 2
    rounding = 1;
end

bytes = get_file_size(filepath);
[hr_bytes, hr_unit] = human_readable_bytes(bytes);
hr_bytes = round(hr_bytes, rounding);

hr_size = [num2str(hr_bytes) hr_unit];
end


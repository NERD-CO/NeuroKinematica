% Function to randomly trim datasets
function rand_trimmed_data = rand_trim_dataset(data, min_size)
if isempty(data) || ~isnumeric(data)
    rand_trimmed_data = [];
elseif length(data) > min_size
    rand_trimmed_data = randsample(data, min_size);
else
    rand_trimmed_data = data;
end
end
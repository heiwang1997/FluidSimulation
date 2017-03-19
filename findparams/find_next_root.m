function [ next_root ] = find_next_root( arr, this_root, end_iter )

for i = (this_root + 1):(end_iter - 1)
    if (arr(i) * arr(i + 1) <= 0)
        next_root = i;
        break;
    end
end

end


%%% copy a 1D matrix along YZ to form a 3D matrix
function Q = repmat_YZ(x, Y_dimension, Z_dimension)
length_ = length(x);
Q = zeros(length_, Y_dimension, Z_dimension);
for i = 1:Y_dimension
    for j = 1:Z_dimension
        Q(:, i, j) = x;
    end
end

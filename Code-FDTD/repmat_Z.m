%%% copy a 2D matrix along Z to form a 3D matrix
function Q=repmat_Z(x,Z_dimension)
size_=size(x);
Q=zeros(size_(1),size_(2),Z_dimension);
for i=1:Z_dimension
    Q(:,:,i)=x;
end

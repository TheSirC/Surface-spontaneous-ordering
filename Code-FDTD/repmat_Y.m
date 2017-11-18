%%% copy a 2D matrix along Y to form a 3D matrix
function Q=repmat_Y(x,Y_dimension)
size_=size(x);
Q=zeros(size_(1),Y_dimension,size_(2));
for i=1:Y_dimension
    Q(:,i,:)=x;
end

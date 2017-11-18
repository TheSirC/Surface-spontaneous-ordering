function Q = Eab2rough2(E_threshold, E_ab, Ez_ab, roughness_function, raduis_Eab2rough2) %#ok<INUSL>
E_ab = 1e-6 * E_ab;
Q = roughness_function;
Q(E_ab >= E_threshold) = 0;
size_roughness_function = size(Q);
center_y = round(size_roughness_function(2)/2);
center_z = round(size_roughness_function(3)/2);
delete_layer = [];
Q_keep = Q;
for i = 1:size_roughness_function(1)
    temp = reshape(Q(i, :, :), size_roughness_function(2), size_roughness_function(3));
    if sum(sum(temp(round(center_y-raduis_Eab2rough2*center_y):round(center_y+raduis_Eab2rough2*center_y), round(center_z-raduis_Eab2rough2*center_z):round(center_z+raduis_Eab2rough2*center_z)))) == 0
        delete_layer = [delete_layer, i]; %#ok<AGROW>
        
    end
end
Q_keep(delete_layer, :, :) = [];
size_Q_keep = size(Q_keep);
Q = ones(size_roughness_function(1), size_roughness_function(2), size_roughness_function(3));
%Q(end:-1:end-size_Q_keep(1)+1,:,:)=Q_keep;
Q(end:-1:end-size_Q_keep(1)+1, :, :) = Q_keep(end:-1:1, :, :);
end
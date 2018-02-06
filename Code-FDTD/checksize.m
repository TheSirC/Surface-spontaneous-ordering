function [ I ] = checksize( image_to_resize, box_size_y, box_size_z, ...
                            conversion_to_real_units_y,conversion_to_real_units_z)
%checksize Function that will ensure that the computation box have at most
%the maximum size for the computation box.
%Otherwise, it will crop the image relative to its center
%A warning will be produced if the latter happens0
% box_size is the maximum size of the computation box in real units

box_size_pix_y = box_size_y/conversion_to_real_units_y;
box_size_pix_z = box_size_z/conversion_to_real_units_z;

if (box_size_pix_y*box_size_pix_z) > 1000
    warning_text = 'The image you have selected is to big\n to be calculated in a resonable amount of time!\n The image will be resized on its center for computation!';
    %Display the warning in the console
    warning(warning_text);
    %and in a dialog box
    f=warndlg(sprintf(warning_text));
    waitfor(f);
    
    %Resize the image
    [Ix, Iy, ~] = size(image_to_resize);
    xymin = min(Ix, Iy);
    dif = (Iy - Ix);
    I = image_to_resize(1:xymin, floor(dif/2):floor(Iy-dif/2)-1);
    
else %Sending back the image if it is not to big
    I = image_to_resize;
end

end
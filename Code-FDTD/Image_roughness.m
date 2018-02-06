function [ roughness ] = Image_roughness( image, roughness_thickness )
%Image_roughness Produces the initial condition for the roughness layer 
%which will be read by the FDTD
%   image is the transformed image retrieved from the MEB
%   rougness_thickness is the amplitude of the 

k = 255/(roughness_thickness-1);
roughness = k*(round(im2uint8(image))./k)+1/2;

end
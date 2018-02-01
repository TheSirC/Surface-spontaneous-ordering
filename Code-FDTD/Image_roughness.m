function [ roughness ] = Image_roughness( image, roughness_thickness )
%Image_roughness Produces the initial condition for the roughness layer 
%which will be read by the FDTD
%   image is the transformed image retrieved from the MEB
%   rougness_thickness is the amplitude of the 

colormap = parula(roughness_thickness);
roughness = dither(image,colormap,roughness_thickness,8);

end
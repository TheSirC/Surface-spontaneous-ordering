# Spontaneous periodic ordering on the surface and in the bulk of dielectrics irradiated by ultrafast laser

## Layout of the project
Entrance point : `FDTD_3D_ripples_LD_main`

## Usage
1. Launch `FDTD_3D_ripples_LD_main`
2. 

## Work to do
- [x] Write the [`Read_plot_image`](https://github.com/TheSirC/Surface-spontaneous-ordering/blob/master/Read_plot_image.m) file to check if able to reconstruct correctly an image without any information on it
- [ ] Develop a function `Image_roughness` which will be the initial condition for the roughness layer which will be read by the FDTD code.
- [ ] Define the amplitude of the roughness layer as a parameter (number of cells) and convert 256 levels of grey into defined number of cells
- [ ] Calculate the field Eab for several images
- [ ] Discuss the location of absorption for several images 

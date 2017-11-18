%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%3D-FDTD code (the function file, the FDTD calucaltion is done here)
%%%to simulate laser pulse propagation in a lossy and
%%%dispersive medium using the Drude model, the transient permittivity change
%%%during the laser pulse is possible to be added in with a dynamically
%%%changing material properties.
%%%Author: Dr.Hao Zhang
%%%Date: March 17, 2014
%%% Lorentz-Drude model added, April 17, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ex,Ey,Ez are in the unit of the FDTD Dennis book, Hx, Hy, Hz are SI units%%%

function [Ez_SI_record_xyzt, Ex_SI_record_xyzt, Ey_SI_record_xyzt, steady_z, steady_y, E_ab, Ex_ab, Ey_ab, Ez_ab, E_0_film, Ex_0_film, Ey_0_film, Ez_0_film, Ex_SI, Ey_SI, Ez_SI, Hx, Hy, Hz, N, epsilon_r_all, air_thickness_act, film_thickness_act, SiO2_thickness_act, ...
    Substract_thickness_act, PML_thickness_act, N_points_x_air, roughness_function, x_corr, y_corr, z_corr, N_points_x_film, N_points_x_SiO2, ...
    N_points_x_Substract, N_points_PML, x_end_film, N_points_x_all, N_points_y_all, N_points_z_all, N_points_y_main, N_points_z_main, E_ab_t, E_Uc_t, E_Ul_t] = ...
    FDTD_3D_ripples_LD_loop ...
    (N_initial, m_, epsilon_r_silicon, epsilon_r_air, epsilon_r_glass, tao_d, f0, omega0, f1, gama1, omega1, f2, gama2, omega2, f3, gama3, omega3, f4, gama4, omega4, ...
    f5, gama5, omega5, order_Snm, air_thickness, film_thickness, SiO2_thickness, ...
    Substract_thickness, PML_thickness, PML_coeff, PML_order, Y_dimension_main, Z_dimension_main, delta_x, laser_calculation_time, T_simulation, add_source, source_edge, points_no_roughness_edge, no_rough_material, ...
    F_av, w0, gau_ord, Tao, lambda, add_roughness_flag, roughness_type, grating_roughness_thick, mean_rough, rough_size, rough_thick, grating_period, grating_width, grating_direction, load_roughness_flag, ...
    roughness_dic, m_1, m_k, Tao_c_l_0, N_cri_phonon, beta, theta_impact, n_2, refreshing_flag, gpu_flag, working_dictionary, ...
    save_dnymics_flag, single_precision_flag, diffusion_flag, pulse_num, Drude_LD_mode, see_depth_animation, see_fields_flag, write_avi_flag, sor_air_edge, record_field_flag, T_av, frequency_field_recording, position_recording)

display('preparing...');
if gpu_flag == 1
    display(sprintf('GPU acceleration enabled'));
    display('reseting gpuDevice...');
    reset(gpuDevice);
end
global COUNTER_T;
%%% some constant %%%
h = 6.626e-34; %%% planck constant J*s%%%
hbar = h / (2 * pi);
eV = 1.602e-19; %%% electron volt J%%%
kB = 1.3807e-23; %%%Boltzmann costant%%% J/K
couloum_e = -1.602e-19; %electric charge of electron
m_e = 9.109e-31; %rest mass of electron kg
epsilon_0 = 8.85e-12; %%% vacuum permittivity F/m
miu_0 = 4 * pi * 1.0e-7; %%% vacuum permeability H/m
c_vacuum = 1 / sqrt(epsilon_0*miu_0); %%% speed of light in vacuum m/s

N_points_x_air = round(air_thickness/delta_x) - 1; %%% the number of cells to be used %%% Because there should be N_cells_x_air points cell, not N_cells_air+1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! But in the code it is N_cells_x_air+1, so we should correct it here.
N_points_x_film = round(film_thickness/delta_x) - 1; %%% Because there should be N_cells_x_film points cell, not N_cells_x_film+1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_points_x_SiO2 = round(SiO2_thickness/delta_x) - 1; %%% Because there should be N_cells_x_SiO2 points cell, not N_cells_x_SiO2+1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_points_x_Substract = round(Substract_thickness/delta_x) - 1; %%% Because there should be N_cells_x_Substract points cell, not N_cells_x_Substract+1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_points_PML = round(PML_thickness/delta_x) - 1; %%% Because there should be N_cells_x_PML points cell, not N_cells_x_PML+1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

air_thickness_act = (N_points_x_air + 1) * delta_x; %%% the actual air thickness
film_thickness_act = (N_points_x_film + 1) * delta_x; %%% the actual film thickness
SiO2_thickness_act = (N_points_x_SiO2 + 1) * delta_x; %%% the actual SiO2 thickness
Substract_thickness_act = (N_points_x_Substract + 1) * delta_x; %%% the actual Substract thickness
PML_thickness_act = (N_points_PML + 1) * delta_x; %%% the actual PML thickness

N_points_y_main = round(Y_dimension_main/delta_x) - 1;

if mod(N_points_y_main+1, 2) == 0
    N_points_y_main = N_points_y_main + 1; %%% make sure there are always odd number of grids on the main grid of Y direction,so that the beam distribution is symetric.
end

Y_dimension_main = (N_points_y_main + 1) * delta_x;

N_points_z_main = round(Z_dimension_main/delta_x) - 1;

if mod(N_points_z_main+1, 2) == 0
    N_points_z_main = N_points_z_main + 1; %%% make sure there are always odd number of grids on the main grid of Y direction,so that the beam distribution is symetric.
end

Z_dimension_main = (N_points_z_main + 1) * delta_x;

N_points_x_all = N_points_x_air + 1 + N_points_x_film + 1 + N_points_x_SiO2 + 1 + N_points_x_Substract + 1 + 2 * (N_points_PML + 1);
N_points_y_all = N_points_y_main + 1 + 2 * (N_points_PML + 1); %%% always an odd integer
N_points_z_all = N_points_z_main + 1 + 2 * (N_points_PML + 1); %%% always an odd integer

f = c_vacuum / lambda; %%% laser frequency Hz %%%
T_cycle = 1 / f; %%% the length of one optical cycle

omega = 2 * pi * f;
k_0 = omega / c_vacuum; %%% wave vector m^-1;

delta_t = (0.5 / (c_vacuum)) * delta_x; %%% time increment
factor_0_5 = c_vacuum * (delta_t / delta_x);

if frequency_field_recording == Inf
    delta_t_recording = delta_t;
else
    delta_t_recording = (T_cycle) / frequency_field_recording; %%% every 1/frequency_field_recording of the optical cycle, recording the E and H field
    %delta_t_recording=(T_cycle)/5;%%% every 1/10 of the optical cycle, recording the E and H field
end

x_start_PML_low = 1;
x_end_PML_low = x_start_PML_low + N_points_PML;
x_start_Substract = x_end_PML_low + 1;
x_end_Substract = x_start_Substract + N_points_x_Substract;
x_start_SiO2 = x_end_Substract + 1;
x_end_SiO2 = x_start_SiO2 + N_points_x_SiO2;
x_start_film = x_end_SiO2 + 1;
x_end_film = x_start_film + N_points_x_film;
x_start_air = x_end_film + 1;
x_end_air = x_start_air + N_points_x_air;
x_start_PML_high = x_end_air + 1;
x_end_PML_high = x_start_PML_high + N_points_PML;

x_corr = (1:N_points_x_all) * delta_x; %%% x corrdinates
y_corr = (1:N_points_y_all) * delta_x; %%% y corrdinates
z_corr = (1:N_points_z_all) * delta_x; %%% z corrdinates

display(sprintf('N_points_y_main=%f', N_points_y_main));
display(sprintf('N_points_z_main=%f', N_points_z_main));
display(sprintf('N_points_x_film=%f', N_points_x_film));

save([working_dictionary, '/save_temp/N_points_PML.mat'], 'N_points_PML');
save([working_dictionary, '/save_temp/N_points_x_all.mat'], 'N_points_x_all');
save([working_dictionary, '/save_temp/N_points_x_air.mat'], 'N_points_x_air');
save([working_dictionary, '/save_temp/N_points_x_SiO2.mat'], 'N_points_x_SiO2');
save([working_dictionary, '/save_temp/N_points_x_Substract.mat'], 'N_points_x_Substract');
save([working_dictionary, '/save_temp/N_points_z_all.mat'], 'N_points_z_all');
save([working_dictionary, '/save_temp/N_points_y_all.mat'], 'N_points_y_all');
save([working_dictionary, '/save_temp/N_points_x_film.mat'], 'N_points_x_film');
save([working_dictionary, '/save_temp/N_points_y_main.mat'], 'N_points_y_main');
save([working_dictionary, '/save_temp/N_points_z_main.mat'], 'N_points_z_main');
save([working_dictionary, '/save_temp/x_corr.mat'], 'x_corr');
save([working_dictionary, '/save_temp/y_corr.mat'], 'y_corr');
save([working_dictionary, '/save_temp/z_corr.mat'], 'z_corr');

N = N_initial * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%%initial electron density distribution    carrier number per m^3%%%
%%%%%%%%%%%%%%%%%%%%
%%% add roughness%%%
%%%%%%%%%%%%%%%%%%%%
roughness_function_first_layer = ones(N_points_y_main+1, N_points_z_main+1); %%% is the b(x,y) in PRB Skolski2012
roughness_function_several_layer = ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %% ajout d'une nouvelle fonction
New2 = [];
New1 = [];
New0 = [];
roughness_function = ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
if add_roughness_flag == 1
    if load_roughness_flag == 1
        roughness_function = importdata([roughness_dic, 'roughness_function.mat']);
    elseif roughness_type == 2 && load_roughness_flag == 0
        mean_rough = mean_rough / rough_size^2;
        roughness_function_first_layer = rand(N_points_y_main+1, N_points_z_main+1);
        
        roughness_function_first_layer(roughness_function_first_layer <= (1 - mean_rough)) = 0;
        roughness_function_first_layer(roughness_function_first_layer > (1 - mean_rough)) = 1;
        [pos_0_y, pos_0_z] = find(roughness_function_first_layer <= (1 - mean_rough));
        [pos_1_y, pos_1_z] = find(roughness_function_first_layer > (1 - mean_rough));
        max_pos_1_y = max(pos_1_y);
        max_pos_1_z = max(pos_1_z);
        for i = 1:length(pos_0_y)
            roughness_function_first_layer(pos_0_y(i), pos_0_z(i)) = 0;
        end
        for i = 1:length(pos_1_y)
            if pos_1_y(i) <= max_pos_1_y - (rough_size - 1) && pos_1_z(i) <= max_pos_1_z - (rough_size - 1)
                roughness_function_first_layer(pos_1_y(i):pos_1_y(i)+(rough_size - 1), pos_1_z(i):pos_1_z(i)+(rough_size - 1)) = 1;
            end
        end
        grating_points = [];
        grating_length_points = round(0.5*grating_period/delta_x);
        temp = 1:grating_length_points:N_points_z_main + 1;
        grating_width = grating_width - grating_period / 2;
        grating_width = round(grating_width/delta_x);
        if mod(length(temp), 2) == 0
            for i = 0:length(temp) / 2 - 1
                grating_points = [grating_points, temp(2*i+1) + grating_width:temp(2*i+2)]; %#ok<AGROW>
            end
        else
            for i = 0:(length(temp) - 1) / 2 - 1
                grating_points = [grating_points, temp(2*i+1) + grating_width:temp(2*i+2)]; %#ok<AGROW>
            end
        end
        
        if grating_direction == 1
            roughness_function_first_layer(:, grating_points) = 0;
        else
            roughness_function_first_layer(grating_points, :) = 0;
        end
        
        for i = 1:rough_thick
            roughness_function(end-(i - 1), :, :) = roughness_function_first_layer;
        end
        
        grating_points_filled = setdiff(1:N_points_z_main+1, grating_points);
        
        
        for i = rough_thick:-1:grating_roughness_thick + 1
            roughness_function(end-(i - 1), :, grating_points_filled) = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Gaussian Case      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif roughness_type == 3 && load_roughness_flag == 0 % Gaussian
        % Parameters to define the surface with a Gaussian spectrum
        mean_rough_gauss = mean_rough / (rough_size^2 * 100);
        %N = N_points_y_main+1;
        rL = 6; %4
        corre_length_x = 200; % in nm
        corre_length_y = 40; % in nm
        roughness_function_several_layer = Discretised_roughness(N_points_y_main+1, rL, mean_rough_gauss, corre_length_x, corre_length_y);
        
        for i = 1:length(roughness_function_several_layer)
            for j = 1:length(roughness_function_several_layer)
                if roughness_function_several_layer(i, j) == 2
                    New2(i, j) = 1;
                    New1(i, j) = 1;
                    New0(i, j) = 1;
                elseif roughness_function_several_layer(i, j) == 1
                    New2(i, j) = 0;
                    New1(i, j) = 1;
                    New0(i, j) = 1;
                else
                    New2(i, j) = 0;
                    New1(i, j) = 0;
                    New0(i, j) = 1;
                end
            end
        end
        
        % Define the upper layers of the sample
        %roughness_function(12,:,:)=New2;
        %roughness_function(11,:,:)=New1;
        %roughness_function(10,:,:)=New0;
        roughness_function(end, :, :) = New2;
        roughness_function(end-1, :, :) = New1;
        %roughness_function(end-2,:,:)=New0;
        mean_roughness = mean(mean(mean(roughness_function(end:-1:end-(rough_thick - 1), points_no_roughness_edge+1:end-points_no_roughness_edge-1, points_no_roughness_edge+1:end-points_no_roughness_edge-1))));
        display(sprintf('mean_roughness_check=%f', mean_roughness));
        
    elseif roughness_type == 1 && load_roughness_flag == 0
        mean_rough = mean_rough / rough_size^2;
        roughness_function_first_layer = rand(N_points_y_main+1, N_points_z_main+1);
        [pos_0_y, pos_0_z] = find(roughness_function_first_layer <= (1 - mean_rough));
        [pos_1_y, pos_1_z] = find(roughness_function_first_layer > (1 - mean_rough));
        max_pos_1_y = max(pos_1_y);
        max_pos_1_z = max(pos_1_z);
        for i = 1:length(pos_0_y)
            roughness_function_first_layer(pos_0_y(i), pos_0_z(i)) = 0;
        end
        for i = 1:length(pos_1_y)
            if pos_1_y(i) <= max_pos_1_y - (rough_size - 1) && pos_1_z(i) <= max_pos_1_z - (rough_size - 1)
                roughness_function_first_layer(pos_1_y(i):pos_1_y(i)+(rough_size - 1), pos_1_z(i):pos_1_z(i)+(rough_size - 1)) = 1;
            end
        end
        
        for i = 1:rough_thick
            roughness_function(end-(i - 1), :, :) = roughness_function_first_layer;
        end
        
    end
    if points_no_roughness_edge ~= 0
        roughness_function_first_layer(:, 1:points_no_roughness_edge) = no_rough_material;
        roughness_function_first_layer(:, end-points_no_roughness_edge:end) = no_rough_material;
        roughness_function_first_layer(1:points_no_roughness_edge, :) = no_rough_material;
        roughness_function_first_layer(end-points_no_roughness_edge:end, :) = no_rough_material;
    end
    
    mean_roughness = mean(mean(mean(roughness_function(end:-1:end-(rough_thick - 1), points_no_roughness_edge+1:end-points_no_roughness_edge-1, points_no_roughness_edge+1:end-points_no_roughness_edge-1))));
    display(sprintf('mean_roughness=%f', mean_roughness));
    N = N .* roughness_function;
    
end

rough_height = zeros(N_points_y_main+1, N_points_z_main+1); %%% acting as the function of an AFM
for i = 1:N_points_y_main + 1
    for j = 1:N_points_z_main + 1
        for k = N_points_x_film + 1:-1:1
            if roughness_function(k, i, j) == 1
                rough_height(i, j) = (N_points_x_film + 1 - k) * delta_x;
                %rough_height(i,j)=k*delta_x;
                break;
            end
        end
    end
end
D_focal_rough = (Y_dimension_main - 2 * source_edge * delta_x) / 2;
gau_ord_rough = 5; %%% super gaussian distribution of order gau_ord, for gau_ord=2, it is a gaussian distribution.
w0_rough = D_focal_rough / (4^(1 / gau_ord)); %%% beam waist at the focus
YZ_rough = meshgrid(y_corr(N_points_PML+1+1:N_points_y_all-(N_points_PML + 1)), z_corr(N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)));
y_shift_rough = y_corr(round(length(y_corr)/2));
z_shift_rough = z_corr(round(length(z_corr)/2));
window_function_rough = exp(-(2 * (abs((YZ_rough - y_shift_rough)).^gau_ord_rough + abs((YZ_rough' - z_shift_rough)).^gau_ord_rough))./(w0_rough^gau_ord_rough));
rough_height = rough_height .* window_function_rough; %%% multiply by a window function to increase the quality of the fft map
rough_height_fft = fftshift(fft2(rough_height, max(2^10, N_points_y_main+1-2*source_edge), max(2^10, N_points_z_main+1-2*source_edge)));

if sor_air_edge ~= 0
    N(:, 1:sor_air_edge, :) = 0;
    N(:, end-sor_air_edge:end, :) = 0;
    N(:, :, 1:sor_air_edge) = 0;
    N(:, :, end-sor_air_edge:end) = 0;
end
wp = ((N * couloum_e^2) ./ (m_ * m_e * epsilon_0)).^0.5; %%% initial plasma frequency rad/s

epsilon_normal_silicon = epsilon_r_silicon .* ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
epsilon_normal_silicon = epsilon_normal_silicon .* roughness_function;
temp = epsilon_normal_silicon(end, :, :);
temp(temp == 0) = 1;
epsilon_normal_silicon(end, :, :) = temp;

sigma_FILM_OPA_silicon = imag(epsilon_normal_silicon) * omega * epsilon_0; %%% effective electric conductivity to account for the one photon absorbtion S/m

nonlinear_relative_permittivity = 0;
epsilon_r_FILM_normal_and_kerr = real(epsilon_normal_silicon) + real(nonlinear_relative_permittivity); %%% the real part of epsilon silicon,inculding the normal and kerr part.
if add_roughness_flag == 1
    epsilon_r_FILM_normal_and_kerr(roughness_function == 0) = 1;
end
if sor_air_edge ~= 0
    epsilon_r_FILM_normal_and_kerr(:, 1:sor_air_edge, :) = 1;
    epsilon_r_FILM_normal_and_kerr(:, end-sor_air_edge:end, :) = 1;
    epsilon_r_FILM_normal_and_kerr(:, :, 1:sor_air_edge) = 1;
    epsilon_r_FILM_normal_and_kerr(:, :, end-sor_air_edge:end) = 1;
end

%%% set up the dielectric constant for the PML layers
epsilon_r_all = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
%epsilon_r_all_0_str=ones(N_points_x_all,N_points_y_all,N_points_z_all); %%% epsilon without strucutre

epsilon_r_all(x_start_PML_low:x_end_PML_low, 1:N_points_PML+1, :) = real(epsilon_r_silicon); %%% the left part of PML layer
epsilon_r_all(x_start_Substract:x_end_Substract, 1:N_points_PML+1, :) = real(epsilon_r_silicon); %%% the left part of PML layer
epsilon_r_all(x_start_SiO2:x_end_SiO2, 1:N_points_PML+1, :) = real(epsilon_r_glass); %%% the left part of PML layer
epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
    (reshape(epsilon_r_FILM_normal_and_kerr(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
epsilon_r_all(x_start_air:x_end_air, 1:N_points_PML+1, :) = real(epsilon_r_air); %%% the left part of PML layer
epsilon_r_all(x_start_PML_high:x_end_PML_high, 1:N_points_PML+1, :) = real(epsilon_r_silicon); %%% the left part of PML layer

epsilon_r_all(x_start_PML_low:x_end_PML_low, N_points_y_all-(N_points_PML):N_points_y_all, :) = real(epsilon_r_silicon); %%% the right part of PML layer
epsilon_r_all(x_start_Substract:x_end_Substract, N_points_y_all-(N_points_PML):N_points_y_all, :) = real(epsilon_r_silicon); %%% the right part of PML layer
epsilon_r_all(x_start_SiO2:x_end_SiO2, N_points_y_all-(N_points_PML):N_points_y_all, :) = real(epsilon_r_glass); %%% the right part of PML layer
epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
    repmat_Y(reshape(epsilon_r_FILM_normal_and_kerr(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
epsilon_r_all(x_start_air:x_end_air, N_points_y_all-(N_points_PML):N_points_y_all, :) = real(epsilon_r_air); %%% the right part of PML layer
epsilon_r_all(x_start_PML_high:x_end_PML_high, N_points_y_all-(N_points_PML):N_points_y_all, :) = real(epsilon_r_silicon); %%% the right part of PML layer

epsilon_r_all(N_points_x_all-(N_points_PML):N_points_x_all, :, :) = real(epsilon_r_air); %%% the upper part of PML layer
epsilon_r_all(1:N_points_PML+1, :, :) = real(epsilon_r_silicon); %%% the lower part of PML layer

%%% set up the dielectric constant for the PML layers, back
epsilon_r_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
    (epsilon_r_FILM_normal_and_kerr(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = real(epsilon_r_silicon);
epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = real(epsilon_r_silicon);

epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ(epsilon_r_FILM_normal_and_kerr(:, 1, 1), N_points_PML+1, N_points_PML+1);
epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
    (epsilon_r_FILM_normal_and_kerr(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);

%%% set up the dielectric constant for the PML layers, front
epsilon_r_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
    (epsilon_r_FILM_normal_and_kerr(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = real(epsilon_r_silicon);
epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = real(epsilon_r_silicon);

epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
    (epsilon_r_FILM_normal_and_kerr(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
    (epsilon_r_FILM_normal_and_kerr(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);

%%% set up the dielectric constant for the PML layers, front

%%% set up the dielectric constant for other parts of the structure

epsilon_r_all(x_start_Substract:x_end_Substract, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = real(epsilon_r_silicon); %%% Substract
epsilon_r_all(x_start_SiO2:x_end_SiO2, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = real(epsilon_r_glass); %%% SiO2
epsilon_r_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
    epsilon_r_FILM_normal_and_kerr;%%% film %%% film
epsilon_r_all(x_start_air:x_end_air, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = real(epsilon_r_air); %%% air

if Drude_LD_mode == 1
    sigma_air = imag(epsilon_r_air) * omega * epsilon_0; %%% electric conductivity S/m
    sigma_glass = imag(epsilon_r_glass) * omega * epsilon_0; %%% electric conductivity S/m
    sigma_substract = imag(epsilon_r_silicon) * omega * epsilon_0; %%% electric conductivity S/m
    sigma_FILM_OPA_silicon = imag(epsilon_normal_silicon) * omega * epsilon_0; %%% effective electric conductivity to account for the one photon absorbtion S/m
    sigma_FILM_excited_silicon = sigma_FILM_OPA_silicon + epsilon_0 * (f0^0.5 * wp).^2 * tao_d; %%% electric conductivity of excited silicon S/m, notice the f0 in the LD model
    
    %%% set up the conducticity for the PML layers
    sigma_all = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
    %sigma_all_0_str=zeros(N_points_x_all,N_points_y_all,N_points_z_all); %%% sigma without structure
    
    sigma_all(x_start_PML_low:x_end_PML_low, 1:N_points_PML+1, :) = sigma_substract; %%% the left part of PML layer
    sigma_all(x_start_Substract:x_end_Substract, 1:N_points_PML+1, :) = sigma_substract; %%% the left part of PML layer
    sigma_all(x_start_SiO2:x_end_SiO2, 1:N_points_PML+1, :) = sigma_glass; %%% the left part of PML layer
    sigma_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
        (reshape(sigma_FILM_excited_silicon(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
    sigma_all(x_start_air:x_end_air, 1:N_points_PML+1, :) = sigma_air; %%% the left part of PML layer
    sigma_all(x_start_PML_high:x_end_PML_high, 1:N_points_PML+1, :) = sigma_air; %%% the left part of PML layer
    
    sigma_all(x_start_PML_low:x_end_PML_low, N_points_y_all-(N_points_PML):N_points_y_all, :) = sigma_substract; %%% the right part of PML layer
    sigma_all(x_start_Substract:x_end_Substract, N_points_y_all-(N_points_PML):N_points_y_all, :) = sigma_substract; %%% the right part of PML layer
    sigma_all(x_start_SiO2:x_end_SiO2, N_points_y_all-(N_points_PML):N_points_y_all, :) = sigma_glass; %%% the right part of PML layer
    sigma_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
        repmat_Y(reshape(sigma_FILM_excited_silicon(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
    
    sigma_all(x_start_air:x_end_air, N_points_y_all-(N_points_PML):N_points_y_all, :) = sigma_air; %%% the right part of PML layer
    sigma_all(x_start_PML_high:x_end_PML_high, N_points_y_all-(N_points_PML):N_points_y_all, :) = sigma_air; %%% the right part of PML layer
    
    sigma_all(N_points_x_all-(N_points_PML):N_points_x_all, :, :) = sigma_air; %%% the upper part of PML layer
    sigma_all(1:N_points_PML+1, :, :) = sigma_substract; %%% the lower part of PML layer
    
    %%% set up the conductivity for the PML layers, back
    sigma_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
        (sigma_FILM_excited_silicon(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
    %sigma_all(x_start_film:x_end_film,1:N_points_PML+1,1:N_points_PML+1)=sigma_FILM_OPA_silicon;
    %sigma_all(x_start_film:x_end_film,N_points_y_all-(N_points_PML):N_points_y_all,1:N_points_PML+1)=sigma_FILM_OPA_silicon;
    
    sigma_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
        (sigma_FILM_excited_silicon(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
    sigma_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
        (sigma_FILM_excited_silicon(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
    
    %%% set up the conductivity for the PML layers, front
    sigma_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
        (sigma_FILM_excited_silicon(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
    %sigma_all(x_start_film:x_end_film,1:N_points_PML+1,N_points_z_all-(N_points_PML):N_points_z_all)=sigma_FILM_OPA_silicon;
    %sigma_all(x_start_film:x_end_film,N_points_y_all-(N_points_PML):N_points_y_all,N_points_z_all-(N_points_PML):N_points_z_all)=sigma_FILM_OPA_silicon;
    sigma_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
        (sigma_FILM_excited_silicon(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
    sigma_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
        (sigma_FILM_excited_silicon(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
    
    %%%
    
    %%% set up the conductivity for other parts of the structure
    sigma_all(x_start_Substract:x_end_Substract, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = sigma_substract; %%% Substract
    sigma_all(x_start_SiO2:x_end_SiO2, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = sigma_glass; %%% SiO2
    sigma_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
        sigma_FILM_excited_silicon;%%% film
    sigma_all(x_start_air:x_end_air, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = sigma_air; %%% air
    
    %%% set up the debye term Chi for the PML layers
    Chi_air = 0; %%% dispersion paramenter, no unit
    Chi_glass = 0; %%% dispersion paramenter, no unit
    Chi_substract = 0; %%% dispersion paramenter, no unit
    Chi_film_excited_silicon = -(f0^0.5 * wp).^2 .* tao_d^2; %%% notice the f0 in the LD model
    
    Chi_all = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
    %Chi_all_0_str=zeros(N_points_x_all,N_points_y_all,N_points_z_all); %%% Chi without structure
    
    Chi_all(x_start_PML_low:x_end_PML_low, 1:N_points_PML+1, :) = Chi_substract; %%% the left part of PML layer
    Chi_all(x_start_Substract:x_end_Substract, 1:N_points_PML+1, :) = Chi_substract; %%% the left part of PML layer
    Chi_all(x_start_SiO2:x_end_SiO2, 1:N_points_PML+1, :) = Chi_glass; %%% the left part of PML layer
    Chi_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
        (reshape(Chi_film_excited_silicon(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
    Chi_all(x_start_air:x_end_air, 1:N_points_PML+1, :) = Chi_air; %%% the left part of PML layer
    Chi_all(x_start_PML_high:x_end_PML_high, 1:N_points_PML+1, :) = Chi_air; %%% the left part of PML layer
    
    Chi_all(x_start_PML_low:x_end_PML_low, N_points_y_all-(N_points_PML):N_points_y_all, :) = Chi_substract; %%% the right part of PML layer
    Chi_all(x_start_Substract:x_end_Substract, N_points_y_all-(N_points_PML):N_points_y_all, :) = Chi_substract; %%% the right part of PML layer
    Chi_all(x_start_SiO2:x_end_SiO2, N_points_y_all-(N_points_PML):N_points_y_all, :) = Chi_glass; %%% the right part of PML layer
    Chi_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
        repmat_Y(reshape(Chi_film_excited_silicon(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
    Chi_all(x_start_air:x_end_air, N_points_y_all-(N_points_PML):N_points_y_all, :) = Chi_air; %%% the right part of PML layer
    Chi_all(x_start_PML_high:x_end_PML_high, N_points_y_all-(N_points_PML):N_points_y_all, :) = Chi_air; %%% the right part of PML layer
    
    Chi_all(N_points_x_all-(N_points_PML):N_points_x_all, :, :) = Chi_air; %%% the upper part of PML layer
    Chi_all(1:N_points_PML+1, :, :) = Chi_substract; %%% the lower part of PML layer
    
    %%% set up Chi for the PML layers, back
    Chi_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
        (Chi_film_excited_silicon(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
    Chi_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
        (Chi_film_excited_silicon(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
    Chi_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
        (Chi_film_excited_silicon(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
    
    %%% set up the conductivity for the PML layers, front
    Chi_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
        (Chi_film_excited_silicon(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
    Chi_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
        (Chi_film_excited_silicon(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
    Chi_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
        (Chi_film_excited_silicon(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
    
    %%%
    
    %%% set up the debye term Chi for other parts of the structure
    Chi_all(x_start_Substract:x_end_Substract, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = Chi_substract; %%% Substract
    Chi_all(x_start_SiO2:x_end_SiO2, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = Chi_glass; %%% SiO2
    Chi_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
        Chi_film_excited_silicon;%%% film
    Chi_all(x_start_air:x_end_air, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = Chi_air; %%% air
    
    %%% set up the debye term t0 for the PML layers
    t0_air = tao_d; %%% dispersion paramenter, s
    t0_film_excited_silicon = tao_d * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%% dispersion paramenter, s
    t0_glass = tao_d; %%% dispersion paramenter, s
    t0_substract = tao_d; %%% dispersion paramenter, s
    
    t0_all = tao_d * ones(N_points_x_all, N_points_y_all, N_points_z_all);
    %t0_all_0_str=tao_d*ones(N_points_x_all,N_points_y_all); %%% t0 without structure
    
    t0_all(x_start_PML_low:x_end_PML_low, 1:N_points_PML+1, :) = t0_substract; %%% the left part of PML layer
    t0_all(x_start_Substract:x_end_Substract, 1:N_points_PML+1, :) = t0_substract; %%% the left part of PML layer
    t0_all(x_start_SiO2:x_end_SiO2, 1:N_points_PML+1, :) = t0_glass; %%% the left part of PML layer
    t0_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
        (reshape(t0_film_excited_silicon(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
    t0_all(x_start_air:x_end_air, 1:N_points_PML+1, :) = t0_air; %%% the left part of PML layer
    t0_all(x_start_PML_high:x_end_PML_high, 1:N_points_PML+1, :) = t0_air; %%% the left part of PML layer
    
    t0_all(x_start_PML_low:x_end_PML_low, N_points_y_all-(N_points_PML):N_points_y_all, :) = t0_substract; %%% the right part of PML layer
    t0_all(x_start_Substract:x_end_Substract, N_points_y_all-(N_points_PML):N_points_y_all, :) = t0_substract; %%% the right part of PML layer
    t0_all(x_start_SiO2:x_end_SiO2, N_points_y_all-(N_points_PML):N_points_y_all, :) = t0_glass; %%% the right part of PML layer
    t0_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
        repmat_Y(reshape(t0_film_excited_silicon(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
    t0_all(x_start_air:x_end_air, N_points_y_all-(N_points_PML):N_points_y_all, :) = t0_air; %%% the right part of PML layer
    t0_all(x_start_PML_high:x_end_PML_high, N_points_y_all-(N_points_PML):N_points_y_all, :) = t0_air; %%% the right part of PML layer
    
    t0_all(N_points_x_all-(N_points_PML):N_points_x_all, :, :) = t0_air; %%% the upper part of PML layer
    t0_all(1:N_points_PML+1, :, :) = t0_substract; %%% the lower part of PML layer
    
    %%% set up Chi for the PML layers, back
    t0_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
        (t0_film_excited_silicon(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
    %t0_all(x_start_film:x_end_film,1:N_points_PML+1,1:N_points_PML+1)=tao_d;
    %t0_all(x_start_film:x_end_film,N_points_y_all-(N_points_PML):N_points_y_all,1:N_points_PML+1)=tao_d;
    t0_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
        (t0_film_excited_silicon(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
    t0_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
        (t0_film_excited_silicon(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
    %%% set up the conductivity for the PML layers, front
    t0_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
        (t0_film_excited_silicon(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
    %t0_all(x_start_film:x_end_film,1:N_points_PML+1,N_points_z_all-(N_points_PML):N_points_z_all)=tao_d;
    %t0_all(x_start_film:x_end_film,N_points_y_all-(N_points_PML):N_points_y_all,N_points_z_all-(N_points_PML):N_points_z_all)=tao_d;
    t0_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
        (t0_film_excited_silicon(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
    t0_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
        (t0_film_excited_silicon(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
    %%%
    
    %%% set up the debye term t0 for other parts of the structure
    t0_all(x_start_Substract:x_end_Substract, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = t0_substract; %%% Substract
    t0_all(x_start_SiO2:x_end_SiO2, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = t0_glass; %%% SiO2
    t0_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
        t0_film_excited_silicon;%%% film
    t0_all(x_start_air:x_end_air, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = t0_air; %%% air
elseif Drude_LD_mode == 2
    gama0 = 1 ./ tao_d;
    [S01_all, S01_film_excited_silicon] = Snm(1, f0, gama0, omega0, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S12 for the PML layers
    [S02_all, S02_film_excited_silicon] = Snm(2, f0, gama0, omega0, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S13 for the PML layers
    [S03_all, S03_film_excited_silicon] = Snm(3, f0, gama0, omega0, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
elseif Drude_LD_mode == 3
    gama0 = 1 ./ tao_d;
    [S01_all, S01_film_excited_silicon] = Snm(1, f0, gama0, omega0, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S12 for the PML layers
    [S02_all, S02_film_excited_silicon] = Snm(2, f0, gama0, omega0, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S13 for the PML layers
    [S03_all, S03_film_excited_silicon] = Snm(3, f0, gama0, omega0, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    %%% set up the lorentz-Drude term S11 for the PML layers
    [S11_all, S11_film_excited_silicon] = Snm(1, f1, gama1, omega1, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S12 for the PML layers
    [S12_all, S12_film_excited_silicon] = Snm(2, f1, gama1, omega1, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S13 for the PML layers
    [S13_all, S13_film_excited_silicon] = Snm(3, f1, gama1, omega1, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S21 for the PML layers
    [S21_all, S21_film_excited_silicon] = Snm(1, f2, gama2, omega2, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S22 for the PML layers
    [S22_all, S22_film_excited_silicon] = Snm(2, f2, gama2, omega2, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S23 for the PML layers
    [S23_all, S23_film_excited_silicon] = Snm(3, f2, gama2, omega2, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S21 for the PML layers
    [S31_all, S31_film_excited_silicon] = Snm(1, f3, gama3, omega3, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S22 for the PML layers
    [S32_all, S32_film_excited_silicon] = Snm(2, f3, gama3, omega3, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S23 for the PML layers
    [S33_all, S33_film_excited_silicon] = Snm(3, f3, gama3, omega3, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S21 for the PML layers
    [S41_all, S41_film_excited_silicon] = Snm(1, f4, gama4, omega4, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S22 for the PML layers
    [S42_all, S42_film_excited_silicon] = Snm(2, f4, gama4, omega4, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S23 for the PML layers
    [S43_all, S43_film_excited_silicon] = Snm(3, f4, gama4, omega4, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S21 for the PML layers
    [S51_all, S51_film_excited_silicon] = Snm(1, f5, gama5, omega5, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S22 for the PML layers
    [S52_all, S52_film_excited_silicon] = Snm(2, f5, gama5, omega5, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
    
    %%% set up the lorentz-Drude term S23 for the PML layers
    [S53_all, S53_film_excited_silicon] = Snm(3, f5, gama5, omega5, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
        x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm);
end

n_half_cycle = round(laser_calculation_time/(T_av * T_cycle)); %%% how many half cycle contained in the pulse
laser_calculation_time = n_half_cycle * (T_av * T_cycle);
n_delta_t_half_cycle = (T_av * T_cycle) / delta_t; %%% how many delta_t contained in half a cycle
n_delta_t_av = round(n_delta_t_half_cycle); %%% for calculate the transient I_gen, because E_0^2=2(Integral(E_0^2*cos(omega*t)^2,0,T/2));
delta_t_I_gen = n_delta_t_av * delta_t; %%% the delta_t in refreshing I_gen
zz = delta_t_I_gen / delta_x^2; %%% for the diffusion equation

delta_t_recording_ = round(delta_t_recording/delta_t);
delta_t_recording = delta_t_recording_ * delta_t; %%% every ~ 1/frequency_field_recording of the optical cycle, recording the E field
delta_counter_t_recording = delta_t_recording / delta_t;
t_points_record_all = 0;
t = 0; COUNTER_T = 0; %%% start time
t_steps_total_pulse = 0;
while t < T_simulation
    if t < laser_calculation_time
        t_steps_total_pulse = t_steps_total_pulse + 1;
        if mod(COUNTER_T, delta_counter_t_recording) == 0
            t_points_record_all = t_points_record_all + 1;
        end
    end
    t = t + delta_t; COUNTER_T = COUNTER_T + 1;
end
t_steps_total = COUNTER_T; %%% total number of time steps

if record_field_flag == 1
    Ez_SI_record_xyzt = zeros(N_points_x_all-(2 * (N_points_PML + 1)), N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge), n_delta_t_av, 'single'); %%% record the Ez field
    Ex_SI_record_xyzt = zeros(N_points_x_all-(2 * (N_points_PML + 1)), N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge), n_delta_t_av, 'single'); %%% record the Ez field
    Ey_SI_record_xyzt = zeros(N_points_x_all-(2 * (N_points_PML + 1)), N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge), n_delta_t_av, 'single'); %%% record the Ez field
    %Ez_SI_record_xyt=zeros(N_points_x_all-(2*(N_points_PML+1)),N_points_y_all-(2*(N_points_PML+1)+2*source_edge),t_points_record,'single'); %%% record the Ez field
    %Ez_SI_record_xzt=zeros(N_points_x_all-(2*(N_points_PML+1)),N_points_z_all-(2*(N_points_PML+1)+2*source_edge),t_points_record,'single'); %%% record the Ez field
    
end
t_points_record = 0;
source_pos_x = x_start_air + round((air_thickness)/delta_x);
%source_pos_x=round(N_points_x_all/2);
source_pos_x_center = source_pos_x(1) + (length(source_pos_x) - 1) / 2; %%% the center of the beam
source_pos_y = N_points_PML + 1 + 1 + source_edge:N_points_PML + 1 + N_points_y_main + 1 - source_edge;
source_pos_y_center = source_pos_y(1) + (length(source_pos_y) - 1) / 2; %%% the center of the beam
source_pos_z = N_points_PML + 1 + 1 + source_edge:N_points_PML + 1 + N_points_z_main + 1 - source_edge;
source_pos_z_center = source_pos_z(1) + (length(source_pos_z) - 1) / 2; %%% the center of the beam
length_source_pos_y = length(source_pos_y);
length_source_pos_z = length(source_pos_z);
F_r_cm = zeros(1, length(source_pos_y), length(source_pos_z));
for i = 1:length(source_pos_y)
    for j = 1:length(source_pos_z)
        F_r_cm(1, i, j) = 2 * F_av * exp(-(2 * (abs((source_pos_y(i) - source_pos_y_center)*delta_x).^gau_ord + abs((source_pos_z(j) - source_pos_z_center)*delta_x).^gau_ord))./(w0^gau_ord)); %%% fluence dependence on r %%%
    end
end
I0 = (1e4 * F_r_cm) / ((1 / ((4 * log(2)) / pi)^0.5) * Tao); %%%peak laser intensity %%% J*s^-1*m^-2 calculated by F=Integrate[I0t]dt

if single_precision_flag == 1
    if gpu_flag == 1
        I0 = gpuArray(single(I0));
    else
        I0 = single(I0);
    end
elseif single_precision_flag == 0
    if gpu_flag == 1
        I0 = gpuArray(I0);
    end
end

display(sprintf('delta_x= %fnm', 1e9*delta_x));
display(sprintf('air_thickness_act= %fnm', 1e9*air_thickness_act));
display(sprintf('film_thickness_act= %fnm', 1e9*film_thickness_act));
display(sprintf('SiO2_thickness_act= %fnm', 1e9*SiO2_thickness_act));
display(sprintf('Substract_thickness_act= %fnm', 1e9*Substract_thickness_act));
display(sprintf('Y_dimension_main=%f micron', Y_dimension_main*1e6));

display(sprintf('%d cells air', N_points_x_air+1));
display(sprintf('%d cells film', N_points_x_film+1));
display(sprintf('%d cells SiO2', N_points_x_SiO2+1));
display(sprintf('%d cells Substract', N_points_x_Substract+1));
display(sprintf('%d cells PML', N_points_PML+1));

if pulse_num == 0 || pulse_num == 1
    figure('name', 'epsilon_r_all xy plane')
    pcolor(1e6*y_corr, 1e6*x_corr, epsilon_r_all(:, :, round(N_points_z_all/2))); xlabel('y (\mum)'); ylabel('x (\mum)'); shading interp;
    
    figure('name', 'epsilon_r_all_xz_plane')
    pcolor(1e6*z_corr, 1e6*x_corr, reshape(epsilon_r_all(:, round(N_points_z_all/2), :), N_points_x_all, N_points_z_all)); xlabel('\mum'); ylabel('\mum'); shading interp;
    
    figure('name', 'epsilon_r_all yz plane,surface of the film')
    pcolor(1e6*z_corr, 1e6*y_corr, reshape(epsilon_r_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
    axis equal;
    
    figure('name', 'epsilon_r_all_all yz plane,middle of the film')
    pcolor(1e6*z_corr, 1e6*y_corr, reshape(epsilon_r_all(round((x_end_film + x_start_film)/2), :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
    axis equal;
    
    if Drude_LD_mode == 1
        figure('name', 'sigma_all_xy_plane')
        pcolor(1e6*y_corr, 1e6*x_corr, sigma_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'sigma_all_xz_plane')
        pcolor(1e6*z_corr, 1e6*x_corr, reshape(sigma_all(:, round(N_points_z_all/2), :), N_points_x_all, N_points_z_all)); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'sigma_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(sigma_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
        
        figure('name', 'Chi_all')
        pcolor(1e6*y_corr, 1e6*x_corr, Chi_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'Chi_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(Chi_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
        
        figure('name', 't0_all')
        pcolor(1e6*y_corr, 1e6*x_corr, t0_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
    elseif Drude_LD_mode == 2
        figure('name', 'S01_all_xy_plane')
        pcolor(1e6*y_corr, 1e6*x_corr, S01_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S01_all_xz_plane')
        pcolor(1e6*z_corr, 1e6*x_corr, reshape(S01_all(:, round(N_points_z_all/2), :), N_points_x_all, N_points_z_all)); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S01_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(S01_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
        
        figure('name', 'S02_all')
        pcolor(1e6*y_corr, 1e6*x_corr, S02_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S02_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(S02_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
        
        figure('name', 'S03_all_xy_plane')
        pcolor(1e6*y_corr, 1e6*x_corr, S03_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S03_all_xz_plane')
        pcolor(1e6*z_corr, 1e6*x_corr, reshape(S03_all(:, round(N_points_z_all/2), :), N_points_x_all, N_points_z_all)); xlabel('\mum'); ylabel('\mum'); shading interp;
    elseif Drude_LD_mode == 3
        figure('name', 'S11_all')
        pcolor(1e6*y_corr, 1e6*x_corr, S11_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S11_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(S11_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
        
        figure('name', 'S12_all')
        pcolor(1e6*y_corr, 1e6*x_corr, S12_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S12_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(S12_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
        
        figure('name', 'S13_all_xy_plane')
        pcolor(1e6*y_corr, 1e6*x_corr, S13_all(:, :, round(N_points_z_all/2))); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S13_all_xz_plane')
        pcolor(1e6*z_corr, 1e6*x_corr, reshape(S13_all(:, round(N_points_z_all/2), :), N_points_x_all, N_points_z_all)); xlabel('\mum'); ylabel('\mum'); shading interp;
        
        figure('name', 'S13_all yz plane,surface of the film')
        pcolor(1e6*z_corr, 1e6*y_corr, reshape(S13_all(x_end_film, :, :), N_points_y_all, N_points_z_all)); xlabel('z (\mum)'); ylabel('y (\mum)'); shading interp;
        axis equal;
    end
else
    figure('name', 'roughness_function xy plane'); pause;
    pcolor(roughness_function(:, :, round(N_points_z_main/2))); shading interp;
    
    figure('name', 'roughness_function xz plane'); pause;
    pcolor(reshape(roughness_function(:, round(N_points_y_main/2), :), N_points_x_film+1, N_points_z_main+1)); shading interp;
    
    figure('name', 'roughness_function first layer')
    pcolor(reshape(roughness_function(end, :, :), N_points_y_main+1, N_points_z_main+1)); shading interp
    
    figure('name', 'roughness_function second layer')
    pcolor(reshape(roughness_function(end-1, :, :), N_points_y_main+1, N_points_z_main+1)); shading interp
    
    figure('name', 'roughness_function third layer')
    pcolor(reshape(roughness_function(end-2, :, :), N_points_y_main+1, N_points_z_main+1)); shading interp
    
    figure('name', 'roughness_function fourth layer')
    pcolor(reshape(roughness_function(end-3, :, :), N_points_y_main+1, N_points_z_main+1)); shading interp
    
    figure('name', 'roughness_function fifth layer')
    pcolor(reshape(roughness_function(end-4, :, :), N_points_y_main+1, N_points_z_main+1)); shading interp
    
    figure('name', 'roughness_function sixth layer')
    pcolor(reshape(roughness_function(end-5, :, :), N_points_y_main+1, N_points_z_main+1)); shading interp
    
end
display(sprintf('grid size: %f X %f X %f', N_points_x_all, N_points_y_all, N_points_z_all));

if single_precision_flag == 1
    display(sprintf('single precision enabled'));
end

if Drude_LD_mode == 1
    display(sprintf('using Drude model (sigma,Chi,t0 form)'));
elseif Drude_LD_mode == 2
    display(sprintf('using Drude model (LD form but with omega0=0)'));
elseif Drude_LD_mode == 3
    display(sprintf('using Lorentz-Drude model'));
end
display(sprintf('press any key to start the FDTD main loop...press Ctrl+C to terminate'));
if pulse_num == 0
    pause;
end
display('preparing...');
close all;
pause(0.1);
%gd = gpuDevice();
%%% initialize matrix for the field

if single_precision_flag == 1
    if gpu_flag == 0
        epsilon_r_all = single(epsilon_r_all);
        if Drude_LD_mode == 1
            sigma_all = single(sigma_all);
            Chi_all = single(Chi_all);
            t0_all = single(t0_all);
            I_Ez = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S_Ez = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            I_Ex = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            I_Ey = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S_Ex = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S_Ey = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        elseif Drude_LD_mode == 2
            S01_all = single(S01_all);
            S02_all = single(S02_all);
            S03_all = single(S03_all);
            S0_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        elseif Drude_LD_mode == 3
            S01_all = single(S01_all);
            S02_all = single(S02_all);
            S03_all = single(S03_all);
            S11_all = single(S11_all);
            S12_all = single(S12_all);
            S13_all = single(S13_all);
            S21_all = single(S21_all);
            S22_all = single(S22_all);
            S23_all = single(S23_all);
            S31_all = single(S31_all);
            S32_all = single(S32_all);
            S33_all = single(S33_all);
            S41_all = single(S41_all);
            S42_all = single(S42_all);
            S43_all = single(S43_all);
            S51_all = single(S51_all);
            S52_all = single(S52_all);
            S53_all = single(S53_all);
            S0_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S1_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S1_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S1_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S2_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S2_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S2_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S3_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S3_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S3_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S4_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S4_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S4_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S5_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S5_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S5_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        end
        Dz = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Dz field
        Hx = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Hx field
        Hy = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Hy field
        
        I_Hx = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');
        I_Hy = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');
        I_Hz = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');
        
        Hz = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Hz field
        Dx = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Dx field
        Dy = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Dy field
        
        I_Dx = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');
        I_Dy = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');
        I_Dz = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');
        
        Ex = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Ey = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Ez = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
    elseif gpu_flag == 1 %%% use CUDA
        epsilon_r_all = gpuArray(single(epsilon_r_all));
        if Drude_LD_mode == 1
            sigma_all = gpuArray(single(sigma_all));
            Chi_all = gpuArray(single(Chi_all));
            t0_all = gpuArray(single(t0_all));
            I_Ez = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S_Ez = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            I_Ex = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            I_Ey = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S_Ex = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S_Ey = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        elseif Drude_LD_mode == 2
            S01_all = gpuArray(single(S01_all));
            S02_all = gpuArray(single(S02_all));
            S03_all = gpuArray(single(S03_all));
            S0_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        elseif Drude_LD_mode == 3
            S01_all = gpuArray(single(S01_all));
            S02_all = gpuArray(single(S02_all));
            S03_all = gpuArray(single(S03_all));
            S11_all = gpuArray(single(S11_all));
            S12_all = gpuArray(single(S12_all));
            S13_all = gpuArray(single(S13_all));
            S21_all = gpuArray(single(S21_all));
            S22_all = gpuArray(single(S22_all));
            S23_all = gpuArray(single(S23_all));
            S31_all = gpuArray(single(S31_all));
            S32_all = gpuArray(single(S32_all));
            S33_all = gpuArray(single(S33_all));
            S41_all = gpuArray(single(S41_all));
            S42_all = gpuArray(single(S42_all));
            S43_all = gpuArray(single(S43_all));
            S51_all = gpuArray(single(S51_all));
            S52_all = gpuArray(single(S52_all));
            S53_all = gpuArray(single(S53_all));
            S0_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S0_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S1_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S1_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S1_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S2_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S2_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S2_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S3_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S3_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S3_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S4_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S4_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S4_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S5_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S5_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
            S5_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        end
        Dz = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Dz field
        Hx = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Hx field
        Hy = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Hy field
        
        I_Hx = gpuArray.zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');
        I_Hy = gpuArray.zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');
        I_Hz = gpuArray.zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');
        
        Hz = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Hz field
        Dx = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Dx field
        Dy = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single'); %%% Dy field
        
        I_Dx = gpuArray.zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');
        I_Dy = gpuArray.zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');
        I_Dz = gpuArray.zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');
        
        Ex = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Ey = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Ez = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
    end
elseif single_precision_flag == 0
    if gpu_flag == 0
        if Drude_LD_mode == 1
            I_Ez = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S_Ez = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            I_Ex = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            I_Ey = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S_Ex = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S_Ey = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            
        elseif Drude_LD_mode == 2
            S0_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        elseif Drude_LD_mode == 3
            S0_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S1_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S1_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S1_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S2_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S2_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S2_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S3_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S3_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S3_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S4_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S4_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S4_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S5_x = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S5_y = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S5_z = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        end
        Dz = zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Dz field
        Hx = zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Hx field
        Hy = zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Hy field
        
        I_Hx = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1);
        I_Hy = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1);
        I_Hz = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all);
        
        Hz = zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Hz field
        Dx = zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Dx field
        Dy = zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Dy field
        
        I_Dx = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1);
        I_Dy = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1);
        I_Dz = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all);
        Ex = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Ey = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Ez = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
    elseif gpu_flag == 1 %%% use CUDA
        epsilon_r_all = gpuArray(epsilon_r_all);
        if Drude_LD_mode == 1
            sigma_all = gpuArray(sigma_all);
            Chi_all = gpuArray(Chi_all);
            t0_all = gpuArray(t0_all);
            I_Ez = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S_Ez = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            I_Ex = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            I_Ey = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S_Ex = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S_Ey = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        elseif Drude_LD_mode == 2
            S01_all = gpuArray(S01_all);
            S02_all = gpuArray(S02_all);
            S03_all = gpuArray(S03_all);
            S0_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        elseif Drude_LD_mode == 3
            S01_all = gpuArray(S01_all);
            S02_all = gpuArray(S02_all);
            S03_all = gpuArray(S03_all);
            S11_all = gpuArray(S11_all);
            S12_all = gpuArray(S12_all);
            S13_all = gpuArray(S13_all);
            S21_all = gpuArray(S21_all);
            S22_all = gpuArray(S22_all);
            S23_all = gpuArray(S23_all);
            S31_all = gpuArray(S31_all);
            S32_all = gpuArray(S32_all);
            S33_all = gpuArray(S33_all);
            S41_all = gpuArray(S41_all);
            S42_all = gpuArray(S42_all);
            S43_all = gpuArray(S43_all);
            S51_all = gpuArray(S51_all);
            S52_all = gpuArray(S52_all);
            S53_all = gpuArray(S53_all);
            S0_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S0_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S1_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S1_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S1_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S2_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S2_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S2_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S3_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S3_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S3_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S4_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S4_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S4_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S5_x = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S5_y = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
            S5_z = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        end
        Dz = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Dz field
        Hx = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Hx field
        Hy = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Hy field
        
        I_Hx = gpuArray.zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1);
        I_Hy = gpuArray.zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1);
        I_Hz = gpuArray.zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all);
        
        Hz = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Hz field
        Dx = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Dx field
        Dy = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all); %%% Dy field
        
        I_Dx = gpuArray.zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1);
        I_Dy = gpuArray.zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1);
        I_Dz = gpuArray.zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all);
        
        Ex = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Ey = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Ez = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
    end
end

if gpu_flag == 0
    if single_precision_flag == 1
        Integral_Ex_square = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Integral_Ey_square = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Integral_Ez_square = zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
    elseif single_precision_flag == 0
        Integral_Ex_square = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Integral_Ey_square = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Integral_Ez_square = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
    end
elseif gpu_flag == 1 %%% use CUDA
    if single_precision_flag == 1
        Integral_Ex_square = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Integral_Ey_square = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
        Integral_Ez_square = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all, 'single');
    elseif single_precision_flag == 0
        Integral_Ex_square = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Integral_Ey_square = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
        Integral_Ez_square = gpuArray.zeros(N_points_x_all, N_points_y_all, N_points_z_all);
    end
end

counter_I_gen = 0;

E_ab = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
Ex_ab = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
Ey_ab = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
Ez_ab = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
E_ab_t = zeros(1, n_half_cycle+1); %%% total absorbed laser energy (J/cm^2) at each time
E_Uc_t = zeros(1, n_half_cycle+1); %%% total carrier energy (J/cm^2) at each time
E_Ul_t = zeros(1, n_half_cycle+1); %%% total lattice energy (J/cm^2) at each time

[row_rough_1, col_rough_1] = find(roughness_function_first_layer == 1);
[row_rough_0, col_rough_0] = find(roughness_function_first_layer == 0);

E_0_film = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%% the amplitude of E field in the film V/m
F_gen_film = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%% fluence
if refreshing_flag == 1
    Tc = 300 * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%% initial carrier temperautre K
    m_ = m_1 + m_k * Tc;
    Uc = zeros(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%% initial carrier total energy density J/cm^3
    Tl = 300 * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%%initial lattice temperature K
    Cc = 3 * N * kB; %%% carrier heat capacity %%% J/m^3 K from classical statistical mechanics
    v_c_mean = ((3 * kB * Tc) ./ (m_ * m_e)).^0.5; %%% the mean velocity of carrier, from classical idea gas m/s
    Kc = (1 / 3) * Cc .* v_c_mean.^2 * tao_d; %%% carrier thermal conducticity W/m K
    
    D0 = ((kB * Tc .* tao_d) ./ (m_ * m_e)); %%% carrier diffusivity m^2/s %%% from Einstin relation and the relation between mobility and scattering time
    %D0=0*((kB*Tc*tao_d)./(m_*m_e));%%% carrier diffusivity m^2/s
    Eg = 1.12 * eV * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1); %%% bandgap energy J
    Cl = 1e6 * (1.978 + 3.54e-4 * Tl - 3.68 * Tl.^-2); %%% initial lattice heat capacity J/m^3 K
    Kl = 1e2 * (1585 * Tl.^-1.23); %%% initial lattice thermal conductivity W/m K
    Ul = Cl .* Tl; %%% initial lattice energy density J/m^3%%%
    Tao_c_l = Tao_c_l_0 * (1 + (N / N_cri_phonon).^2);
    alfa_OPA_silicon = (4 * pi * imag(epsilon_r_silicon.^0.5) / lambda); %%% initial linear absorption coefficient at 800nm m^-1
end

n_normal_silicon = epsilon_r_silicon.^0.5;

Gai_3_real = n_2 ./ (3 ./ (4 * real(n_normal_silicon).^2 * epsilon_0 * (c_vacuum))); %%% third order optical susceptibility real part (m^2/V^2)
Gai_3_im = 1i * ((2 * real(n_normal_silicon).^2 * (c_vacuum) * epsilon_0) / (3 * k_0)) * beta; %%% (m^2/V^2) third order optical susceptibility imaginary part, fomule in SI unit
Gai_3 = Gai_3_real + Gai_3_im; %%% third order optical susceptibility

nonlinear_relative_permittivity = (3 / 4) * Gai_3 .* E_0_film.^2;
nonlinear_relative_permittivity = nonlinear_relative_permittivity .* roughness_function;
if Drude_LD_mode == 3
    epsilon_excited_silicon = epsilon_normal_silicon + (f0 * wp.^2) ./ (omega0^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama0) + ...
        (f1 * wp.^2) ./ (omega1^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama1) + ...
        (f2 * wp.^2) ./ (omega2^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama2) + ...
        (f3 * wp.^2) ./ (omega3^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama3) + ...
        (f4 * wp.^2) ./ (omega4^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama4) + ...
        (f5 * wp.^2) ./ (omega5^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama5) + ...
        +nonlinear_relative_permittivity;%%% delectric constant of excited silicon
elseif Drude_LD_mode == 1
    epsilon_excited_silicon = epsilon_normal_silicon - ((wp / (2 * pi * f)).^2) .* (1 ./ (1 + 1i * (1 ./ (2 * pi * f * tao_d)))) + nonlinear_relative_permittivity; %%% delectric constant of excited silicon
elseif Drude_LD_mode == 2
    epsilon_excited_silicon = epsilon_normal_silicon + (f0 * wp.^2) ./ (omega0^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama0) + nonlinear_relative_permittivity; %%% delectric constant of excited silicon
end
n_excited_silicon = epsilon_excited_silicon.^0.5;
if refreshing_flag == 0
    display(sprintf('epsion_exci=%f+%fi', real(epsilon_excited_silicon(1, sor_air_edge+1, sor_air_edge+1)), imag(epsilon_excited_silicon(1, sor_air_edge+1, sor_air_edge+1))));
    wavelentgh_div_delta_x = ((lambda / real(n_excited_silicon(1, sor_air_edge+1, sor_air_edge+1)))) / delta_x;
    display(sprintf('wavelentgh_div_delta_x=%f', wavelentgh_div_delta_x));
end
alfa_excited_silicon = (4 * pi * imag(n_excited_silicon) / lambda); %%% initial linear absorption coefficient at 800nm m^-1

%%% intitialize matrix for the PML layer
xn_temp = PML_coeff * ((1:N_points_PML + 1) / (N_points_PML + 1)).^PML_order;
xn = [turn_around(xn_temp), zeros(1, N_points_x_all-2*(N_points_PML + 1)), xn_temp];
fi1 = xn;
gi1 = xn;
gi2 = 1 ./ (1 + xn);
gi3 = (1 - xn) ./ (1 + xn);
fi2 = gi2;
fi3 = gi3;

xn_temp1 = 1 * PML_coeff * ((1:N_points_PML + 1) / (N_points_PML + 1)).^PML_order;
xn_temp2 = 1 * PML_coeff * ((1:N_points_PML) / (N_points_PML)).^PML_order;
xn_left = [turn_around(xn_temp1), zeros(1, N_points_x_all-2*(N_points_PML + 1)), xn_temp2]; %%% crucial important

fi1_left = xn_left;
gi2_left = 1 ./ (1 + xn_left); %%% crucial important
gi3_left = (1 - xn_left) ./ (1 + xn_left); %%% crucial important
fi2_left = gi2_left;
fi3_left = gi3_left;

xn_right = [turn_around(xn_temp2), zeros(1, N_points_x_all-2*(N_points_PML + 1)), xn_temp1]; %%% crucial important
fi1_right = xn_right;
gi2_right = 1 ./ (1 + xn_right); %%% crucial important
gi3_right = (1 - xn_right) ./ (1 + xn_right); %%% crucial important
fi2_right = gi2_right;
fi3_right = gi3_right;

yn_temp = PML_coeff * ((1:N_points_PML + 1) / (N_points_PML + 1)).^PML_order;
yn = [turn_around(yn_temp), zeros(1, N_points_y_all-2*(N_points_PML + 1)), yn_temp];
fj1 = yn;
gj1 = yn;
gj2 = 1 ./ (1 + yn);
gj3 = (1 - yn) ./ (1 + yn);
fj2 = gj2;
fj3 = gj3;

yn_temp1 = 1 * PML_coeff * ((1:N_points_PML + 1) / (N_points_PML + 1)).^PML_order;
yn_temp2 = 1 * PML_coeff * ((1:N_points_PML) / (N_points_PML)).^PML_order;
yn_left = [turn_around(yn_temp1), zeros(1, N_points_y_all-2*(N_points_PML + 1)), yn_temp2]; %%% crucial important

fj1_left = yn_left;
gj2_left = 1 ./ (1 + yn_left); %%% crucial important
gj3_left = (1 - yn_left) ./ (1 + yn_left); %%% crucial important
fj2_left = gj2_left;
fj3_left = gj3_left;

yn_right = [turn_around(yn_temp2), zeros(1, N_points_y_all-2*(N_points_PML + 1)), yn_temp1]; %%% crucial important
fj1_right = yn_right;
gj2_right = 1 ./ (1 + yn_right); %%% crucial important
gj3_right = (1 - yn_right) ./ (1 + yn_right); %%% crucial important
fj2_right = gj2_right;
fj3_right = gj3_right;

zn_temp = PML_coeff * ((1:N_points_PML + 1) / (N_points_PML + 1)).^PML_order;
zn = [turn_around(zn_temp), zeros(1, N_points_z_all-2*(N_points_PML + 1)), zn_temp];
fk1 = zn;
gk1 = zn;
gk2 = 1 ./ (1 + zn);
gk3 = (1 - zn) ./ (1 + zn);
fk2 = gk2;
fk3 = gk3;

zn_temp1 = 1 * PML_coeff * ((1:N_points_PML + 1) / (N_points_PML + 1)).^PML_order;
zn_temp2 = 1 * PML_coeff * ((1:N_points_PML) / (N_points_PML)).^PML_order;
zn_left = [turn_around(zn_temp1), zeros(1, N_points_z_all-2*(N_points_PML + 1)), zn_temp2]; %%% crucial important

fk1_left = zn_left;
gk2_left = 1 ./ (1 + zn_left); %%% crucial important
gk3_left = (1 - zn_left) ./ (1 + zn_left); %%% crucial important
fk2_left = gk2_left;
fk3_left = gk3_left;

zn_right = [turn_around(zn_temp2), zeros(1, N_points_z_all-2*(N_points_PML + 1)), zn_temp1]; %%% crucial important
fk1_right = zn_right;
gk2_right = 1 ./ (1 + zn_right); %%% crucial important
gk3_right = (1 - zn_right) ./ (1 + zn_right); %%% crucial important
fk2_right = gk2_right;
fk3_right = gk3_right;

gj3gk3_3D = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');
gj2gk2_3D = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');
gi1_3D = zeros(N_points_x_all, N_points_y_all-1, N_points_z_all-1, 'single');

for i = 1:N_points_x_all
    for j = 1:N_points_y_all - 1
        for k = 1:N_points_z_all - 1
            gj3gk3_3D(i, j, k) = gj3_left(j) * gk3_right(k);
            gj2gk2_3D(i, j, k) = gj2_left(j) * gk2_right(k);
            gi1_3D(i, j, k) = gi1(i);
        end
    end
end

if gpu_flag == 1
    gj3gk3_3D = gpuArray(gj3gk3_3D);
    gj2gk2_3D = gpuArray(gj2gk2_3D);
    gi1_3D = gpuArray(gi1_3D);
end

gi3gk3_3D = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');
gi2gk2_3D = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');
gj1_3D = zeros(N_points_x_all-1, N_points_y_all, N_points_z_all-1, 'single');

for i = 1:N_points_x_all - 1
    for j = 1:N_points_y_all
        for k = 1:N_points_z_all - 1
            gi3gk3_3D(i, j, k) = gi3_left(i) * gk3_right(k);
            gi2gk2_3D(i, j, k) = gi2_left(i) * gk2_right(k);
            gj1_3D(i, j, k) = gj1(j);
        end
    end
end

if gpu_flag == 1
    gi3gk3_3D = gpuArray(gi3gk3_3D);
    gi2gk2_3D = gpuArray(gi2gk2_3D);
    gj1_3D = gpuArray(gj1_3D);
end

gi3gj3_3D = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');
gi2gj2_3D = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');
gk1_3D = zeros(N_points_x_all-1, N_points_y_all-1, N_points_z_all, 'single');

for i = 1:N_points_x_all - 1
    for j = 1:N_points_y_all - 1
        for k = 1:N_points_z_all
            gi3gj3_3D(i, j, k) = gi3_left(i) * gj3_right(j);
            gi2gj2_3D(i, j, k) = gi2_left(i) * gj2_right(j);
            gk1_3D(i, j, k) = gk1(k);
        end
    end
end

if gpu_flag == 1
    gi3gj3_3D = gpuArray(gi3gj3_3D);
    gi2gj2_3D = gpuArray(gi2gj2_3D);
    gk1_3D = gpuArray(gk1_3D);
end

gj2gk2_3D = factor_0_5 * gj2gk2_3D; %%% multiply factor_0_5 here other than in the while loop to increase speed a little little bit
gi2gk2_3D = factor_0_5 * gi2gk2_3D; %%% multiply factor_0_5 here other than in the while loop to increase speed a little little bit
gi2gj2_3D = factor_0_5 * gi2gj2_3D; %%% multiply factor_0_5 here other than in the while loop to increase speed a little little bit

steady_z = [];
steady_y = [];
t = 0; COUNTER_T = 0;
if Drude_LD_mode == 1
    exp_delta_t_dvi_t0_all = exp(-(delta_t)./(t0_all));
    delta_t_dvi_to_all = (delta_t ./ t0_all);
    delta_t_div_eps_0 = delta_t / epsilon_0;
end
right_x = int16(2:N_points_x_all);
right_y = int16(2:N_points_y_all);
right_z = int16(2:N_points_z_all);
left_x = int16(1:N_points_x_all-1);
left_y = int16(1:N_points_y_all-1);
left_z = int16(1:N_points_z_all-1);
Subs1 = substruct('()', {':', right_y, right_z});
Subs2 = substruct('()', {':', left_y, right_z});
Subs3 = substruct('()', {':', right_y, left_z});
Subs4 = substruct('()', {':', left_y, left_z});

Subs5 = substruct('()', {right_x, ':', right_z});
Subs6 = substruct('()', {right_x, ':', left_z});
Subs7 = substruct('()', {left_x, ':', right_z});
Subs8 = substruct('()', {left_x, ':', left_z});

Subs9 = substruct('()', {right_x, right_y, ':'});
Subs10 = substruct('()', {right_x, left_y, ':'});
Subs11 = substruct('()', {left_x, right_y, ':'});
Subs12 = substruct('()', {left_x, left_y, ':'});

pre_S0_x = 0;
pre_S0_y = 0;
pre_S0_z = 0;
pre_S1_x = 0;
pre_S1_y = 0;
pre_S1_z = 0;
pre_S2_x = 0;
pre_S2_y = 0;
pre_S2_z = 0;
pre_S3_x = 0;
pre_S3_y = 0;
pre_S3_z = 0;
pre_S4_x = 0;
pre_S4_y = 0;
pre_S4_z = 0;
pre_S5_x = 0;
pre_S5_y = 0;
pre_S5_z = 0;

if see_fields_flag == 1
    h_animation = figure('name', 'animation');
    if write_avi_flag == 1
        avi_obj = avifile([working_dictionary, 'animation.avi']); %#ok<DAVIFL>
    end
end

display('start FDTD loop');
while t < T_simulation
    tic;
    %Calculate the Dx field
    %Calculate the Dx field
    %Calculate the Dx field
    
    %%matrix form
    curl = Hz(:, 2:N_points_y_all, 2:N_points_z_all) - Hz(:, 1:N_points_y_all-1, 2:N_points_z_all) - Hy(:, 2:N_points_y_all, 2:N_points_z_all) + Hy(:, 2:N_points_y_all, 1:N_points_z_all-1);
    %curl=subsref(Hz-Hy,Subs1)-subsref(Hz,Subs2)+subsref(Hy,Subs3); %%% explicitly calling subsref, a little bit faster, but still much slower than not doing subscript indexing, the build in function subsref may has many overheads.
    I_Dx = I_Dx + curl;
    Dx(:, 2:N_points_y_all, 2:N_points_z_all) = gj3gk3_3D .* Dx(:, 2:N_points_y_all, 2:N_points_z_all) + ...
        gj2gk2_3D .* (curl + gi1_3D .* I_Dx);
    
    %Calculate the Dy field
    %Calculate the Dy field
    %Calculate the Dy field
    
    %%matrix form
    curl = Hx(2:N_points_x_all, :, 2:N_points_z_all) - Hx(2:N_points_x_all, :, 1:N_points_z_all-1) - Hz(2:N_points_x_all, :, 2:N_points_z_all) + Hz(1:N_points_x_all-1, :, 2:N_points_z_all);
    %curl=subsref(Hx-Hz,Subs5)-subsref(Hx,Subs6)+subsref(Hz,Subs7); %%% explicitly calling subsref, a little bit faster, but still much slower than not doing subscript indexing, the build in function subsref may has many overheads.
    I_Dy = I_Dy + curl;
    Dy(2:N_points_x_all, :, 2:N_points_z_all) = gi3gk3_3D .* Dy(2:N_points_x_all, :, 2:N_points_z_all) + ...
        gi2gk2_3D .* (curl + gj1_3D .* I_Dy);
    
    %Calculate the Dz field
    %Calculate the Dz field
    %Calculate the Dz field
    
    %%matrix form
    curl = Hy(2:N_points_x_all, 2:N_points_y_all, :) - Hy(1:N_points_x_all-1, 2:N_points_y_all, :) - Hx(2:N_points_x_all, 2:N_points_y_all, :) + Hx(2:N_points_x_all, 1:N_points_y_all-1, :);
    %curl=subsref(Hy-Hx,Subs9)-subsref(Hy,Subs11)+subsref(Hx,Subs10); %%% explicitly calling subsref, a little bit faster, but still much slower than not doing subscript indexing, the build in function subsref may has many overheads.
    I_Dz = I_Dz + curl;
    Dz(2:N_points_x_all, 2:N_points_y_all, :) = gi3gj3_3D .* Dz(2:N_points_x_all, 2:N_points_y_all, :) + ...
        gi2gj2_3D .* (curl + gk1_3D .* I_Dz);
    if add_source == 1
        if t < laser_calculation_time
            % Add source
            I0t = I0 .* exp(-4*log(2)*((t - (0.5 * laser_calculation_time)) / Tao)^2); %%% in W/m^2
            
            %I0t=I0;
            
            %A=((2*I0t)./((reshape(epsilon_r_all(source_pos_x,source_pos_y,source_pos_z),length(source_pos_y),length(source_pos_z))).^0.5*c_vacuum*epsilon_0)).^0.5;%%% incident Ex amplitude,in V/m
            
            A = ((2 * I0t) ./ (1 * c_vacuum * epsilon_0)).^0.5; %%% incident Ex amplitude,in V/m, if the source is in the air
            A = (epsilon_0 / miu_0)^0.5 * A; %%% in the unit of the FDTD book of Dennis
            source = A * cos(omega*t); %%% the wave front at the focal plane is flat
            
            %source=(epsilon_0/miu_0)^0.5*((2*(I0.*exp(-4*log(2)*((t-(0.5*laser_calculation_time))/Tao)^2)))./(1*c_vacuum*epsilon_0)).^0.5*cos(omega*t);
            
            %Dz(source_pos_x,source_pos_y,source_pos_z)=Dz(source_pos_x,source_pos_y,source_pos_z)+epsilon_r_all(source_pos_x,source_pos_y,source_pos_z)...
            %    .*reshape(source,1,length_source_pos_y,length_source_pos_z);
            
            %Dy(source_pos_x,source_pos_y,source_pos_z)=Dy(source_pos_x,source_pos_y,source_pos_z)+epsilon_r_all(source_pos_x,source_pos_y,source_pos_z)...
            %   .*reshape(source,1,length_source_pos_y,length_source_pos_z);
            
            Dz(source_pos_x, source_pos_y, source_pos_z) = Dz(source_pos_x, source_pos_y, source_pos_z) + source; %%% if the source is in the air
            %Dy(source_pos_x,source_pos_y,source_pos_z)=Dy(source_pos_x,source_pos_y,source_pos_z)+source; %%% if the source is in the air
        end
    end
    if Drude_LD_mode == 3
        %%% calculate Ex field
        %%% calculate Ex field
        %%% calculate Ex field
        
        pre_pre_S0_x = pre_S0_x;
        pre_S0_x = S0_x;
        pre_pre_S1_x = pre_S1_x;
        pre_S1_x = S1_x;
        pre_pre_S2_x = pre_S2_x;
        pre_S2_x = S2_x;
        pre_pre_S3_x = pre_S3_x;
        pre_S3_x = S3_x;
        pre_pre_S4_x = pre_S4_x;
        pre_S4_x = S4_x;
        pre_pre_S5_x = pre_S5_x;
        pre_S5_x = S5_x;
        S0_x = S01_all .* S0_x - S02_all .* pre_pre_S0_x + S03_all .* Ex;
        S1_x = S11_all .* S1_x - S12_all .* pre_pre_S1_x + S13_all .* Ex;
        S2_x = S21_all .* S2_x - S22_all .* pre_pre_S2_x + S23_all .* Ex;
        S3_x = S31_all .* S3_x - S32_all .* pre_pre_S3_x + S33_all .* Ex;
        S4_x = S41_all .* S4_x - S42_all .* pre_pre_S4_x + S43_all .* Ex;
        S5_x = S51_all .* S5_x - S52_all .* pre_pre_S5_x + S53_all .* Ex;
        Ex = (Dx - S0_x - S1_x - S2_x - S3_x - S4_x - S5_x) ./ epsilon_r_all;
        
        %%% calculate Ey field
        %%% calculate Ey field
        %%% calculate Ey field
        
        pre_pre_S0_y = pre_S0_y;
        pre_S0_y = S0_y;
        pre_pre_S1_y = pre_S1_y;
        pre_S1_y = S1_y;
        pre_pre_S2_y = pre_S2_y;
        pre_S2_y = S2_y;
        pre_pre_S3_y = pre_S3_y;
        pre_S3_y = S3_y;
        pre_pre_S4_y = pre_S4_y;
        pre_S4_y = S4_y;
        pre_pre_S5_y = pre_S5_y;
        pre_S5_y = S5_y;
        S0_y = S01_all .* S0_y - S02_all .* pre_pre_S0_y + S03_all .* Ey;
        S1_y = S11_all .* S1_y - S12_all .* pre_pre_S1_y + S13_all .* Ey;
        S2_y = S21_all .* S2_y - S22_all .* pre_pre_S2_y + S23_all .* Ey;
        S3_y = S31_all .* S3_y - S32_all .* pre_pre_S3_y + S33_all .* Ey;
        S4_y = S41_all .* S4_y - S42_all .* pre_pre_S4_y + S43_all .* Ey;
        S5_y = S51_all .* S5_y - S52_all .* pre_pre_S5_y + S53_all .* Ey;
        Ey = (Dy - S0_y - S1_y - S2_y - S3_y - S4_y - S5_y) ./ epsilon_r_all;
        
        %%% calculate Ez field
        %%% calculate Ez field
        %%% calculate Ez field
        
        pre_pre_S0_z = pre_S0_z;
        pre_S0_z = S0_z;
        pre_pre_S1_z = pre_S1_z;
        pre_S1_z = S1_z;
        pre_pre_S2_z = pre_S2_z;
        pre_S2_z = S2_z;
        pre_pre_S3_z = pre_S3_z;
        pre_S3_z = S3_z;
        pre_pre_S4_z = pre_S4_z;
        pre_S4_z = S4_z;
        pre_pre_S5_z = pre_S5_z;
        pre_S5_z = S5_z;
        S0_z = S01_all .* S0_z - S02_all .* pre_pre_S0_z + S03_all .* Ez;
        S1_z = S11_all .* S1_z - S12_all .* pre_pre_S1_z + S13_all .* Ez;
        S2_z = S21_all .* S2_z - S22_all .* pre_pre_S2_z + S23_all .* Ez;
        S3_z = S31_all .* S3_z - S32_all .* pre_pre_S3_z + S33_all .* Ez;
        S4_z = S41_all .* S4_z - S42_all .* pre_pre_S4_z + S43_all .* Ez;
        S5_z = S51_all .* S5_z - S52_all .* pre_pre_S5_z + S53_all .* Ez;
        Ez = (Dz - S0_z - S1_z - S2_z - S3_z - S4_z - S5_z) ./ epsilon_r_all;
    elseif Drude_LD_mode == 1
        %%% calculate Ex field
        %%% calculate Ex field
        %%% calculate Ex field
        
        %Ex=Dx./epsilon_r_all;
        
        Ex = (Dx - I_Ex - exp_delta_t_dvi_t0_all .* S_Ex) ./ (epsilon_r_all + (sigma_all * delta_t_div_eps_0) + Chi_all .* delta_t_dvi_to_all);
        I_Ex = I_Ex + (sigma_all * delta_t_div_eps_0) .* Ex;
        S_Ex = exp_delta_t_dvi_t0_all .* S_Ex + Chi_all .* delta_t_dvi_to_all .* Ex;
        
        %Ey=Dy./epsilon_r_all;
        
        Ey = (Dy - I_Ey - exp_delta_t_dvi_t0_all .* S_Ey) ./ (epsilon_r_all + (sigma_all * delta_t_div_eps_0) + Chi_all .* delta_t_dvi_to_all);
        I_Ey = I_Ey + (sigma_all * delta_t_div_eps_0) .* Ey;
        S_Ey = exp_delta_t_dvi_t0_all .* S_Ey + Chi_all .* delta_t_dvi_to_all .* Ey;
        
        %%% calculate Ez field
        %%% calculate Ez field
        %%% calculate Ez field
        
        %Ez=Dz./epsilon_r_all;
        
        Ez = (Dz - I_Ez - exp_delta_t_dvi_t0_all .* S_Ez) ./ (epsilon_r_all + (sigma_all * delta_t_div_eps_0) + Chi_all .* delta_t_dvi_to_all);
        I_Ez = I_Ez + (sigma_all * delta_t_div_eps_0) .* Ez;
        S_Ez = exp_delta_t_dvi_t0_all .* S_Ez + Chi_all .* delta_t_dvi_to_all .* Ez;
    elseif Drude_LD_mode == 2
        %%% calculate Ex field
        %%% calculate Ex field
        %%% calculate Ex field
        
        %Ex=Dx./epsilon_r_all;
        
        pre_pre_S0_x = pre_S0_x;
        pre_S0_x = S0_x;
        S0_x = S01_all .* S0_x - S02_all .* pre_pre_S0_x + S03_all .* Ex;
        Ex = (Dx - S0_x) ./ epsilon_r_all;
        
        %Ey=Dy./epsilon_r_all;
        
        pre_pre_S0_y = pre_S0_y;
        pre_S0_y = S0_y;
        S0_y = S01_all .* S0_y - S02_all .* pre_pre_S0_y + S03_all .* Ey;
        Ey = (Dy - S0_y) ./ epsilon_r_all;
        
        %%% calculate Ez field
        %%% calculate Ez field
        %%% calculate Ez field
        
        %Ez=Dz./epsilon_r_all;
        
        pre_pre_S0_z = pre_S0_z;
        pre_S0_z = S0_z;
        S0_z = S01_all .* S0_z - S02_all .* pre_pre_S0_z + S03_all .* Ez;
        Ez = (Dz - S0_z) ./ epsilon_r_all;
    end
    %%% Remember: part of the PML is E=0 at the edges, PEC
    %%% Remember: part of the PML is E=0 at the edges, PEC
    %%% Remember: part of the PML is E=0 at the edges, PEC
    
    Ex(1, :, :) = 0;
    Ey(1, :, :) = 0;
    Ez(1, :, :) = 0;
    
    Ex(N_points_x_all, :, :) = 0;
    Ey(N_points_x_all, :, :) = 0;
    Ez(N_points_x_all, :, :) = 0;
    
    Ex(:, 1, :) = 0;
    Ey(:, 1, :) = 0;
    Ez(:, 1, :) = 0;
    
    Ex(:, N_points_y_all, :) = 0;
    Ey(:, N_points_y_all, :) = 0;
    Ez(:, N_points_y_all, :) = 0;
    
    Ex(:, :, 1) = 0;
    Ey(:, :, 1) = 0;
    Ez(:, :, 1) = 0;
    
    Ex(:, :, N_points_z_all) = 0;
    Ey(:, :, N_points_z_all) = 0;
    Ez(:, :, N_points_z_all) = 0;
    
    %Calculate the Hx field
    %Calculate the Hx field
    %Calculate the Hx field
    
    %%matrix form
    curl = Ey(:, 1:N_points_y_all-1, 2:N_points_z_all) - Ey(:, 1:N_points_y_all-1, 1:N_points_z_all-1) - Ez(:, 2:N_points_y_all, 1:N_points_z_all-1) + Ez(:, 1:N_points_y_all-1, 1:N_points_z_all-1);
    %curl=subsref(Ey,Subs2)-subsref(Ez,Subs3)+subsref(Ez-Ey,Subs4); %%% explicitly calling subsref, a little bit faster, but still much slower than not doing subscript indexing, the build in function subsref may has many overheads.
    I_Hx = I_Hx + curl;
    Hx(:, 1:N_points_y_all-1, 1:N_points_z_all-1) = gj3gk3_3D .* Hx(:, 1:N_points_y_all-1, 1:N_points_z_all-1) + ...
        gj2gk2_3D .* (curl + gi1_3D .* I_Hx);
    
    %Calculate the Hy field
    %Calculate the Hy field
    %Calculate the Hy field
    
    %%% matrix form
    curl = Ez(2:N_points_x_all, :, 1:N_points_z_all-1) - Ez(1:N_points_x_all-1, :, 1:N_points_z_all-1) - Ex(1:N_points_x_all-1, :, 2:N_points_z_all) + Ex(1:N_points_x_all-1, :, 1:N_points_z_all-1);
    %curl=subsref(Ez,Subs6)-subsref(Ex,Subs7)+subsref(Ex-Ez,Subs8); %%% explicitly calling subsref, a little bit faster, but still much slower than not doing subscript indexing, the build in function subsref may has many overheads.
    I_Hy = I_Hy + curl;
    Hy(1:N_points_x_all-1, :, 1:N_points_z_all-1) = gi3gk3_3D .* Hy(1:N_points_x_all-1, :, 1:N_points_z_all-1) + ...
        gi2gk2_3D .* (curl + gj1_3D .* I_Hy);
    
    %Calculate the Hz field
    %Calculate the Hz field
    %Calculate the Hz field
    
    %%% matrix form
    curl = Ex(1:N_points_x_all-1, 2:N_points_y_all, :) - Ex(1:N_points_x_all-1, 1:N_points_y_all-1, :) - Ey(2:N_points_x_all, 1:N_points_y_all-1, :) + Ey(1:N_points_x_all-1, 1:N_points_y_all-1, :);
    %curl=subsref(Ex,Subs11)-subsref(Ey,Subs10)+subsref(Ey-Ex,Subs12); %%% explicitly calling subsref, a little bit faster, but still much slower than not doing subscript indexing, the build in function subsref may has many overheads.
    I_Hz = I_Hz + curl;
    Hz(1:N_points_x_all-1, 1:N_points_y_all-1, :) = gi3gj3_3D .* Hz(1:N_points_x_all-1, 1:N_points_y_all-1, :) + ...
        gi2gj2_3D .* (curl + gk1_3D .* I_Hz);
    
    Ex_SI = Ex * (miu_0 / epsilon_0)^0.5; %%% in V/m
    Ey_SI = Ey * (miu_0 / epsilon_0)^0.5; %%% in V/m
    Ez_SI = Ez * (miu_0 / epsilon_0)^0.5; %%% in V/m
    
    Integral_Ex_square = Integral_Ex_square + Ex_SI.^2 * delta_t;
    Integral_Ey_square = Integral_Ey_square + Ey_SI.^2 * delta_t;
    Integral_Ez_square = Integral_Ez_square + Ez_SI.^2 * delta_t;
    
    if gpu_flag == 1
        wait(gd);
    end
    time_used_inside = toc();
    
    if mod(COUNTER_T, n_delta_t_av) == 0
        display('updating field amplitude and draw the animation...');
        counter_I_gen = counter_I_gen + 1;
        Ex_0_square = Integral_Ex_square / ((T_cycle / 2) * T_av); %%% amplitude square V^2/m^2
        Ey_0_square = Integral_Ey_square / ((T_cycle / 2) * T_av); %%% amplitude square V^2/m^2
        Ez_0_square = Integral_Ez_square / ((T_cycle / 2) * T_av); %%% amplitude square V^2/m^2
        E_0_square = Ex_0_square + Ey_0_square + Ez_0_square;
        
        %E_0_square_film=E_0_square(x_start_film:x_end_film,(N_points_PML+1+1):N_points_y_all-(N_points_PML+1),N_points_PML+1+1:N_points_z_all-(N_points_PML+1)); %%% amplitude square V^2/m^2
        E_0_film = E_0_square(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)).^0.5; %%% amplitude V/m
        Ex_0_film = Ex_0_square(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)).^0.5; %%% amplitude V/m
        Ey_0_film = Ey_0_square(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)).^0.5; %%% amplitude V/m
        Ez_0_film = Ez_0_square(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)).^0.5; %%% amplitude V/m
        
        I_gen_film = double(0.5*c_vacuum.*real(n_excited_silicon).*epsilon_0.*(E_0_film).^2); %%% optical intensity in the film W/m^2
        Q_energy = 0.5 * c_vacuum * (epsilon_0) * alfa_excited_silicon .* real(n_excited_silicon) .* E_0_film.^2; %%% energy dissipated per second in (J/m^3)*s^-1
        Q_energy_x = 0.5 * c_vacuum * (epsilon_0) * alfa_excited_silicon .* real(n_excited_silicon) .* Ex_0_film.^2;
        Q_energy_y = 0.5 * c_vacuum * (epsilon_0) * alfa_excited_silicon .* real(n_excited_silicon) .* Ey_0_film.^2;
        Q_energy_z = 0.5 * c_vacuum * (epsilon_0) * alfa_excited_silicon .* real(n_excited_silicon) .* Ez_0_film.^2;
        F_gen_film = F_gen_film + I_gen_film .* delta_t_I_gen; %%% fluence in the film J/m^2
        
        %figure(9);mesh(reshape(I_gen_film(end,:,:),N_points_y_all-2*(N_points_PML+1),N_points_z_all-2*(N_points_PML+1)));
        Integral_Ex_square = 0 * Integral_Ex_square;
        Integral_Ey_square = 0 * Integral_Ey_square;
        Integral_Ez_square = 0 * Integral_Ez_square;
        
        if refreshing_flag == 1
            %%% solving carrier generation and diffusion equation %%%
            %%% solving carrier generation and diffusion equation %%%
            %%% solving carrier generation and diffusion equation %%%
            if diffusion_flag == 1
                %%% for the filled part %%%
                for ii_1 = 1:length(row_rough_1)
                    i = row_rough_1(ii_1);
                    j = col_rough_1(ii_1);
                    %%% generating tridiagonal matrix for solving a linear system of equation for x direction
                    a1 = [zz * D0(1, i, j)', (zz / 4) * (D0(3:N_points_x_film+1, i, j)' + 4 * D0(2:N_points_x_film, i, j)' - D0(1:N_points_x_film-1, i, j)'), zz * D0(N_points_x_film+1, i, j)']; %%% matrix elements
                    b1 = 2 * zz * D0(:, i, j)'; %%% matrix elements
                    c1 = [-zz * D0(1, i, j)', (zz / 4) * (D0(3:N_points_x_film+1, i, j)' - 4 * D0(2:N_points_x_film, i, j)' - D0(1:N_points_x_film-1, i, j)'), -zz * D0(N_points_x_film+1, i, j)']; %%% matrix elements
                    diagonal1_x = 1 + b1;
                    low1_x = [c1(2:N_points_x_film), c1(N_points_x_film+1) - a1(N_points_x_film-1)];
                    up1_x = [c1(1) - a1(1), -a1(2:N_points_x_film)];
                    %factor_N_rate=(2e29-N)./2e29;
                    %f_x_t_1=factor_N_rate(:,i,j)'.*((alfa_OPA_silicon.*I_gen_film(:,i,j)')/(h*f)+(beta*(I_gen_film(:,i,j)').^2)/(2*h*f)+theta_impact.*I_gen_film(:,i,j)'.*N(:,i,j)');%%% consider impact ionization
                    f_x_t_1 = ((alfa_OPA_silicon .* I_gen_film(:, i, j)') / (h * f) + (beta * (I_gen_film(:, i, j)').^2) / (2 * h * f) + theta_impact .* I_gen_film(:, i, j)' .* N(:, i, j)'); %%% consider impact ionization
                    N(:, i, j) = (tri_diag_2(up1_x, diagonal1_x, low1_x, N(:, i, j)'+delta_t_I_gen*f_x_t_1))';
                end
            else
                for ii_1 = 1:length(row_rough_1)
                    i = row_rough_1(ii_1);
                    j = col_rough_1(ii_1);
                    %%% generating tridiagonal matrix for solving a linear system of equation for x direction
                    %factor_N_rate=(2e29-N)./2e29;
                    %f_x_t_1=factor_N_rate(:,i,j)'.*((alfa_OPA_silicon.*I_gen_film(:,i,j)')/(h*f)+(beta*(I_gen_film(:,i,j)').^2)/(2*h*f)+theta_impact.*I_gen_film(:,i,j)'.*N(:,i,j)');%%% consider impact ionization
                    f_x_t_1 = ((alfa_OPA_silicon .* I_gen_film(:, i, j)) / (h * f) + (beta * (I_gen_film(:, i, j)).^2) / (2 * h * f) + theta_impact .* I_gen_film(:, i, j) .* N(:, i, j)); %%% consider impact ionization
                    %N(:,i,j)=(diag_solver(diagonal1_x,N(:,i,j)+delta_t_I_gen*f_x_t_1));
                    N(:, i, j) = (N(:, i, j) + delta_t_I_gen * f_x_t_1);
                end
                %f_x_t_1=((alfa_OPA_silicon.*I_gen_film(:,row_rough_1,col_rough_1))/(h*f)+(beta*(I_gen_film(:,row_rough_1,col_rough_1)).^2)/(2*h*f)+theta_impact.*I_gen_film(:,row_rough_1,col_rough_1).*N(:,row_rough_1,col_rough_1));%%% consider impact ionization
                %N(:,row_rough_1,col_rough_1)=(N(:,row_rough_1,col_rough_1)+delta_t_I_gen*f_x_t_1);
            end
            %%% for the unfilled part
            if diffusion_flag == 1
                for ii_0 = 1:length(row_rough_0)
                    i = row_rough_0(ii_0);
                    j = col_rough_0(ii_0);
                    a1 = [zz * D0(1, i, j)', (zz / 4) * (D0(3:N_points_x_film+1-1, i, j)' + 4 * D0(2:N_points_x_film-1, i, j)' - D0(1:N_points_x_film-1-1, i, j)'), zz * D0(N_points_x_film+1-1, i, j)']; %%% matrix elements
                    b1 = 2 * zz * D0(1:N_points_x_film, i, j)'; %%% matrix elements
                    c1 = [-zz * D0(1, i, j)', (zz / 4) * (D0(3:N_points_x_film+1-1, i, j)' - 4 * D0(2:N_points_x_film-1, i, j)' - D0(1:N_points_x_film-1-1, i, j)'), -zz * D0(N_points_x_film+1-1, i, j)']; %%% matrix elements
                    diagonal1_x = 1 + b1;
                    low1_x = [c1(2:N_points_x_film-1), c1(N_points_x_film+1-1) - a1(N_points_x_film-1-1)];
                    up1_x = [c1(1) - a1(1), -a1(2:N_points_x_film-1)];
                    %factor_N_rate=(2e29-N)./2e29;
                    %f_x_t_1=factor_N_rate(1:N_points_x_film,i,j)'.*((alfa_OPA_silicon.*I_gen_film(1:N_points_x_film,i,j)')/(h*f)+(beta*(I_gen_film(1:N_points_x_film,i,j)').^2)/(2*h*f)+theta_impact.*I_gen_film(1:N_points_x_film,i,j)'.*N(1:N_points_x_film,i,j)');%%% consider impact ionization
                    f_x_t_1 = ((alfa_OPA_silicon .* I_gen_film(1:N_points_x_film, i, j)') / (h * f) + (beta * (I_gen_film(1:N_points_x_film, i, j)').^2) / (2 * h * f) + theta_impact .* I_gen_film(1:N_points_x_film, i, j)' .* N(1:N_points_x_film, i, j)'); %%% consider impact ionization
                    N(1:N_points_x_film, i, j) = (tri_diag_2(up1_x, diagonal1_x, low1_x, N(1:N_points_x_film, i, j)'+delta_t_I_gen*f_x_t_1))';
                end
            else
                %diagonal1_x=ones(N_points_x_film,1);
                %spmd
                for ii_0 = 1:length(row_rough_0)
                    i = row_rough_0(ii_0);
                    j = col_rough_0(ii_0);
                    %factor_N_rate=(2e29-N)./2e29;
                    %f_x_t_1=factor_N_rate(1:N_points_x_film,i,j)'.*((alfa_OPA_silicon.*I_gen_film(1:N_points_x_film,i,j)')/(h*f)+(beta*(I_gen_film(1:N_points_x_film,i,j)').^2)/(2*h*f)+theta_impact.*I_gen_film(1:N_points_x_film,i,j)'.*N(1:N_points_x_film,i,j)');%%% consider impact ionization
                    f_x_t_1 = ((alfa_OPA_silicon .* I_gen_film(1:N_points_x_film, i, j)) / (h * f) + (beta * (I_gen_film(1:N_points_x_film, i, j)).^2) / (2 * h * f) + theta_impact .* I_gen_film(1:N_points_x_film, i, j) .* N(1:N_points_x_film, i, j)); %%% consider impact ionization
                    %N(1:N_points_x_film,i,j)=(diag_solver(diagonal1_x,N(1:N_points_x_film,i,j)+delta_t_I_gen*f_x_t_1));
                    N(1:N_points_x_film, i, j) = (N(1:N_points_x_film, i, j) + delta_t_I_gen * f_x_t_1);
                end
                %end
                %N=gather(N);
            end
            N = N .* roughness_function;
            Eg = 1.86e-19 - 1.123e-22 * (Tl.^2 ./ (Tl + 1108)) - 2.4e-29 * N.^(1 / 3);
            Eg(Eg < 0) = 0;
            Cc = 3 * N * kB; %%% A N dependent Cc must be updated here, otherwise it valides conservation of energy
            %figure(11);pcolor(reshape(1e-6*N(N_points_x_film+1-1,:,:),N_points_y_main+1,N_points_z_main+1));colormap(jet(2^10));shading interp;
            display(sprintf('N_max=%e cm^{-3}', 1e-6*max(max(max(N)))));
            
            
            %%% solving electron temperature equation %%%
            %%% solving electron temperature equation %%%
            %%% solving electron temperature equation %%%
            
            %%% for the filled part %%%
            if diffusion_flag == 1
                for ii_1 = 1:length(row_rough_1)
                    i = row_rough_1(ii_1);
                    j = col_rough_1(ii_1);
                    a2 = [zz * Kc(1, i, j)', (zz / 4) * (Kc(3:N_points_x_film+1, i, j)' + 4 * Kc(2:N_points_x_film, i, j)' - Kc(1:N_points_x_film-1, i, j)'), zz * Kc(N_points_x_film+1, i, j)'];
                    b2 = 2 * zz * Kc(:, i, j)';
                    c2 = [-zz * Kc(1, i, j)', (zz / 4) * (Kc(3:N_points_x_film+1, i, j)' - 4 * Kc(2:N_points_x_film, i, j)' - Kc(1:N_points_x_film-1, i, j)'), -zz * Kc(N_points_x_film+1, i, j)'];
                    diagonal2_x = Cc(:, i, j)' + b2;
                    low2_x = [c2(2:N_points_x_film), c2(N_points_x_film+1) - a2(N_points_x_film-1)];
                    up2_x = [c2(1) - a2(1), -a2(2:N_points_x_film)];
                    f_x_t_2 = (alfa_excited_silicon(:, i, j)') .* I_gen_film(:, i, j)' - (Cc(:, i, j)' ./ Tao_c_l(:, i, j)') .* (Tc(:, i, j)' - Tl(:, i, j)');
                    Tc(:, i, j) = real((tri_diag_2(up2_x, diagonal2_x, low2_x, (Uc(:, i, j)' - N(:, i, j)' .* Eg(:, i, j)')+delta_t_I_gen*f_x_t_2))');
                end
            else
                for ii_1 = 1:length(row_rough_1)
                    i = row_rough_1(ii_1);
                    j = col_rough_1(ii_1);
                    diagonal2_x = Cc(:, i, j);
                    f_x_t_2 = (alfa_excited_silicon(:, i, j)) .* I_gen_film(:, i, j) - (Cc(:, i, j) ./ Tao_c_l(:, i, j)) .* (Tc(:, i, j) - Tl(:, i, j));
                    %Tc(:,i,j)=(diag_solver(diagonal2_x,(Uc(:,i,j)-N(:,i,j).*Eg(:,i,j))+delta_t_I_gen*f_x_t_2)); %%% diffusion ignored
                    Tc(:, i, j) = ((Uc(:, i, j) - N(:, i, j) .* Eg(:, i, j)) + delta_t_I_gen * f_x_t_2) ./ diagonal2_x;
                end
            end
            
            %%% for the unfilled part %%%
            if diffusion_flag == 1
                for ii_0 = 1:length(row_rough_0)
                    i = row_rough_0(ii_0);
                    j = col_rough_0(ii_0);
                    a2 = [zz * Kc(1, i, j)', (zz / 4) * (Kc(3:N_points_x_film+1-1, i, j)' + 4 * Kc(2:N_points_x_film-1, i, j)' - Kc(1:N_points_x_film-1-1, i, j)'), zz * Kc(N_points_x_film+1-1, i, j)'];
                    b2 = 2 * zz * Kc(1:N_points_x_film, i, j)';
                    c2 = [-zz * Kc(1, i, j)', (zz / 4) * (Kc(3:N_points_x_film+1-1, i, j)' - 4 * Kc(2:N_points_x_film-1, i, j)' - Kc(1:N_points_x_film-1-1, i, j)'), -zz * Kc(N_points_x_film+1-1, i, j)'];
                    diagonal2_x = Cc(1:N_points_x_film, i, j)' + b2;
                    low2_x = [c2(2:N_points_x_film-1), c2(N_points_x_film+1-1) - a2(N_points_x_film-1-1)];
                    up2_x = [c2(1) - a2(1), -a2(2:N_points_x_film-1)];
                    f_x_t_2 = (alfa_excited_silicon(1:N_points_x_film, i, j)') .* I_gen_film(1:N_points_x_film, i, j)' - (Cc(1:N_points_x_film, i, j)' ./ Tao_c_l(1:N_points_x_film, i, j)') .* (Tc(1:N_points_x_film, i, j)' - Tl(1:N_points_x_film, i, j)');
                    Tc(1:N_points_x_film, i, j) = real((tri_diag_2(up2_x, diagonal2_x, low2_x, (Uc(1:N_points_x_film, i, j)' - N(1:N_points_x_film, i, j)' .* Eg(1:N_points_x_film, i, j)')+delta_t_I_gen*f_x_t_2))');
                end
            else
                for ii_0 = 1:length(row_rough_0)
                    i = row_rough_0(ii_0);
                    j = col_rough_0(ii_0);
                    diagonal2_x = Cc(1:N_points_x_film, i, j);
                    f_x_t_2 = (alfa_excited_silicon(1:N_points_x_film, i, j)) .* I_gen_film(1:N_points_x_film, i, j) - (Cc(1:N_points_x_film, i, j) ./ Tao_c_l(1:N_points_x_film, i, j)) .* (Tc(1:N_points_x_film, i, j) - Tl(1:N_points_x_film, i, j));
                    %Tc(1:N_points_x_film,i,j)=(diag_solver(diagonal2_x,(Uc(1:N_points_x_film,i,j)-N(1:N_points_x_film,i,j).*Eg(1:N_points_x_film,i,j))+delta_t_I_gen*f_x_t_2));
                    Tc(1:N_points_x_film, i, j) = ((Uc(1:N_points_x_film, i, j) - N(1:N_points_x_film, i, j) .* Eg(1:N_points_x_film, i, j)) + delta_t_I_gen * f_x_t_2) ./ diagonal2_x;
                end
            end
            Tc(Tc < 300) = 300;
            Tc = Tc .* roughness_function;
            Cc = 3 * N * kB; %%% the refreshing of Cc must before the update of Uc, otherwise it is not stable!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Uc = N .* Eg + Cc .* Tc; %%% carrier energy density J/cm^3
            %figure(10);pcolor(reshape(Tc(N_points_x_film+1-1,:,:),N_points_y_main+1,N_points_z_main+1));zlabel('Tc (K)');colormap(jet(2^10));shading interp;
            display(sprintf('Tc_max=%e K', max(max(max(Tc)))));
            m_ = m_1 + m_k * Tc;
            v_c_mean = ((3 * kB * Tc) ./ (m_ * m_e)).^0.5; %%% the mean velocity of carrier, from classical idea gas m/s
            Kc = (1 / 3) * Cc .* v_c_mean.^2 .* tao_d; %%% carrier thermal conducticity W/m K
            D0 = ((kB * Tc .* tao_d) ./ (m_ * m_e)); %%% carrier diffusivity m^2/s %%% from Einstin relation and the relation between mobility and scattering time
            wp = ((N * couloum_e^2) ./ (m_ * m_e * epsilon_0)).^0.5; %%% initial plasma frequency rad/s
            nonlinear_relative_permittivity = (3 / 4) * Gai_3 .* E_0_film.^2;
            nonlinear_relative_permittivity = nonlinear_relative_permittivity .* roughness_function;
            epsilon_excited_silicon = epsilon_normal_silicon - ((wp / (2 * pi * f)).^2) .* (1 ./ (1 + 1i * (1 ./ (2 * pi * f * tao_d)))) + nonlinear_relative_permittivity; %%% delectric constant of excited silicon
            epsilon_excited_silicon = epsilon_excited_silicon .* roughness_function;
            epsilon_excited_silicon(epsilon_excited_silicon == 0) = 1;
            n_excited_silicon = epsilon_excited_silicon.^0.5;
            alfa_excited_silicon = (4 * pi * imag(n_excited_silicon) / lambda); %%% initial linear absorption coefficient at 800nm m^-1
            Q_energy = 0.5 * c_vacuum * (epsilon_0) * alfa_excited_silicon .* real(n_excited_silicon) .* E_0_film.^2; %%% energy dissipated per second in (J/cm^3)*s^-1
            %%% solving lattice temperature equation %%%
            %%% solving lattice temperature equation %%%
            %%% solving lattice temperature equation %%%
            
            %%% for the filled part %%%
            if diffusion_flag == 1
                for ii_1 = 1:length(row_rough_1)
                    i = row_rough_1(ii_1);
                    j = col_rough_1(ii_1);
                    a3 = [zz * Kl(1, i, j)', (zz / 4) * (Kl(3:N_points_x_film+1, i, j)' + 4 * Kl(2:N_points_x_film, i, j)' - Kl(1:N_points_x_film-1, i, j)'), zz * Kl(N_points_x_film+1, i, j)'];
                    b3 = 2 * zz * Kl(:, i, j)';
                    c3 = [-zz * Kl(1, i, j)', (zz / 4) * (Kl(3:N_points_x_film+1, i, j)' - 4 * Kl(2:N_points_x_film, i, j)' - Kl(1:N_points_x_film-1, i, j)'), -zz * Kl(N_points_x_film+1, i, j)'];
                    diagonal3_x = Cl(:, i, j)' + b3;
                    low3_x = [c3(2:N_points_x_film), c3(N_points_x_film+1) - a3(N_points_x_film-1)];
                    up3_x = [c3(1) - a3(1), -a3(2:N_points_x_film)];
                    f_x_t_3 = (Cc(:, i, j)' ./ Tao_c_l(:, i, j)') .* (Tc(:, i, j)' - Tl(:, i, j)');
                    Tl(:, i, j) = (tri_diag_2(up3_x, diagonal3_x, low3_x, Ul(:, i, j)'+delta_t_I_gen*f_x_t_3))';
                end
            else
                for ii_1 = 1:length(row_rough_1)
                    i = row_rough_1(ii_1);
                    j = col_rough_1(ii_1);
                    diagonal3_x = Cl(:, i, j);
                    f_x_t_3 = (Cc(:, i, j) ./ Tao_c_l(:, i, j)) .* (Tc(:, i, j) - Tl(:, i, j));
                    %Tl(:,i,j)=(diag_solver(diagonal3_x,Ul(:,i,j)+delta_t_I_gen*f_x_t_3));
                    Tl(:, i, j) = (Ul(:, i, j) + delta_t_I_gen * f_x_t_3) ./ diagonal3_x;
                end
            end
            
            %%% for the unfilled part %%%
            if diffusion_flag == 1
                for ii_0 = 1:length(row_rough_0)
                    i = row_rough_0(ii_0);
                    j = col_rough_0(ii_0);
                    a3 = [zz * Kl(1, i, j)', (zz / 4) * (Kl(3:N_points_x_film+1-1, i, j)' + 4 * Kl(2:N_points_x_film-1, i, j)' - Kl(1:N_points_x_film-1-1, i, j)'), zz * Kl(N_points_x_film+1-1, i, j)'];
                    b3 = 2 * zz * Kl(1:N_points_x_film, i, j)';
                    c3 = [-zz * Kl(1, i, j)', (zz / 4) * (Kl(3:N_points_x_film+1-1, i, j)' - 4 * Kl(2:N_points_x_film-1, i, j)' - Kl(1:N_points_x_film-1-1, i, j)'), -zz * Kl(N_points_x_film+1-1, i, j)'];
                    diagonal3_x = Cl(1:N_points_x_film, i, j)' + b3;
                    low3_x = [c3(2:N_points_x_film-1), c3(N_points_x_film+1-1) - a3(N_points_x_film-1-1)];
                    up3_x = [c3(1) - a3(1), -a3(2:N_points_x_film-1)];
                    f_x_t_3 = (Cc(1:N_points_x_film, i, j)' ./ Tao_c_l(1:N_points_x_film, i, j)') .* (Tc(1:N_points_x_film, i, j)' - Tl(1:N_points_x_film, i, j)');
                    Tl(1:N_points_x_film, i, j) = (tri_diag_2(up3_x, diagonal3_x, low3_x, Ul(1:N_points_x_film, i, j)'+delta_t_I_gen*f_x_t_3))';
                end
            else
                for ii_0 = 1:length(row_rough_0)
                    i = row_rough_0(ii_0);
                    j = col_rough_0(ii_0);
                    diagonal3_x = Cl(1:N_points_x_film, i, j);
                    f_x_t_3 = (Cc(1:N_points_x_film, i, j) ./ Tao_c_l(1:N_points_x_film, i, j)) .* (Tc(1:N_points_x_film, i, j) - Tl(1:N_points_x_film, i, j));
                    %Tl(1:N_points_x_film,i,j)=(diag_solver(diagonal3_x,Ul(1:N_points_x_film,i,j)+delta_t_I_gen*f_x_t_3));
                    Tl(1:N_points_x_film, i, j) = (Ul(1:N_points_x_film, i, j) + delta_t_I_gen * f_x_t_3) ./ diagonal3_x;
                end
            end
            
            Tl = Tl .* roughness_function;
            
            Ul = Cl .* Tl; %%% lattice energy density J/cm^3
            %Ucl=Uc+Ul;%%% internal erngy density
            Eg = 1.86e-19 - 1.123e-22 * (Tl.^2 ./ (Tl + 1108)) - 2.4e-29 * N.^(1 / 3);
            Eg(Eg < 0) = 0;
            Tl_for_Cl = Tl;
            %Tl_for_Cl(Tl_for_Cl>1687)=1687;
            Cl = 1e6 * (1.978 + 3.54e-4 * Tl_for_Cl - 3.68 * Tl_for_Cl.^-2); %%% initial lattice heat capacity J/m^3 K
            
            Kl = 1e2 * (1585 * Tl.^-1.23); %%% initial lattice thermal conductivity W/m K
            display(sprintf('Tl_max=%e K', max(max(max(Tl)))));
            
            %%% refreshing parametres %%%
            
            epsilon_r_FILM_normal_and_kerr = real(epsilon_normal_silicon) + real(nonlinear_relative_permittivity); %%% the real part of epsilon silicon,inculding the normal and kerr part.
            sigma_FILM_TPA_silicon = imag(nonlinear_relative_permittivity) * omega * epsilon_0; %%% effective electric conductivity to account for the one photon absorbtion S/m
            
            epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
                (reshape(epsilon_r_FILM_normal_and_kerr(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
            epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
                repmat_Y(reshape(epsilon_r_FILM_normal_and_kerr(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
            epsilon_r_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
                (epsilon_r_FILM_normal_and_kerr(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
            epsilon_r_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
                (epsilon_r_FILM_normal_and_kerr(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
            
            epsilon_r_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
                epsilon_r_FILM_normal_and_kerr;%%% film
            
            sigma_FILM_excited_silicon = sigma_FILM_OPA_silicon + sigma_FILM_TPA_silicon + epsilon_0 * (f0^0.5 * wp).^2 .* tao_d; %%% electric conductivity of excited silicon S/m
            sigma_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
                (reshape(sigma_FILM_excited_silicon(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
            sigma_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
                repmat_Y(reshape(sigma_FILM_excited_silicon(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
            sigma_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
                (sigma_FILM_excited_silicon(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
            sigma_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
                (sigma_FILM_excited_silicon(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
            sigma_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
                sigma_FILM_excited_silicon;%%% film
            
            Chi_film_excited_silicon = -(f0^0.5 * wp).^2 .* tao_d.^2;
            Chi_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
                (reshape(Chi_film_excited_silicon(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
            Chi_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
                repmat_Y(reshape(Chi_film_excited_silicon(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
            Chi_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
                (Chi_film_excited_silicon(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
            Chi_all(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
                (Chi_film_excited_silicon(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
            Chi_all(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
                Chi_film_excited_silicon;%%% film
            
            %%%%%%%%%
            epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
                (epsilon_r_FILM_normal_and_kerr(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
            epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
                (epsilon_r_FILM_normal_and_kerr(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
            epsilon_r_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
                (epsilon_r_FILM_normal_and_kerr(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
            epsilon_r_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
                (epsilon_r_FILM_normal_and_kerr(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
            
            sigma_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
                (sigma_FILM_excited_silicon(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
            sigma_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
                (sigma_FILM_excited_silicon(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
            sigma_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
                (sigma_FILM_excited_silicon(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
            sigma_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
                (sigma_FILM_excited_silicon(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
            
            Chi_all(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
                (Chi_film_excited_silicon(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
            Chi_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
                (Chi_film_excited_silicon(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
            Chi_all(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
                (Chi_film_excited_silicon(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
            Chi_all(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
                (Chi_film_excited_silicon(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
            
        end
        
        temp_steady_z = Ez_ab(end-1, :, :);
        temp_steady_y = Ey_ab(end-1, :, :);
        
        E_ab = E_ab + Q_energy * delta_t_I_gen; %%% absorbed energy density in J/m^3, total
        Ex_ab = Ex_ab + Q_energy_x * delta_t_I_gen; %%% absorbed energy density in J/m^3, x
        Ey_ab = Ey_ab + Q_energy_y * delta_t_I_gen; %%% absorbed energy density in J/m^3, y
        Ez_ab = Ez_ab + Q_energy_z * delta_t_I_gen; %%% absorbed energy density in J/m^3, z
        if refreshing_flag == 1
            E_ab_t(counter_I_gen) = sum(sum(sum(E_ab))) * delta_x^3; %%% absorbed energy density integrated over x and y J
            E_Uc_t(counter_I_gen) = sum(sum(sum(Uc))) * delta_x^3; %%% carrier energy (J)
            E_Ul_t(counter_I_gen) = sum(sum(sum(Ul))) * delta_x^3; %%% lattice energy (J)
        end
        
        N_fft_x = max(2^10, N_points_y_main+1-2*source_edge);
        N_fft_y = max(2^10, N_points_z_main+1-2*source_edge);
        k_0 = 2 * pi / lambda;
        F_spatial = 2 * pi * (1 / (delta_x)) / k_0; %%% sampling frequency (spatial) normliazed to k_0
        F_x = linspace(-F_spatial/2, F_spatial/2, N_fft_x); %%% coordinate in k space %%%
        F_y = linspace(-F_spatial/2, F_spatial/2, N_fft_y); %%% coordinate in k space %%%
        F_xy_mesh = meshgrid(F_x, F_y);
        %%% use to defind the colormap for the fft image
        F_colormap = 1; %%% k_0
        colormap_F = find(abs(F_xy_mesh) > F_colormap);
        color_factor = 1;
        %%% use to defind the colormap for the fft image
        
        E_ab_slice = reshape(E_ab(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        Ex_ab_slice = reshape(Ex_ab(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        Ey_ab_slice = reshape(Ey_ab(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        Ez_ab_slice = reshape(Ez_ab(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        
        E_0_film_slice = reshape(E_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        Ex_0_film_slice = reshape(Ex_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        Ey_0_film_slice = reshape(Ey_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        Ez_0_film_slice = reshape(Ez_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
        
        E_ab_slice_fft = fftshift(fft2(E_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        Ex_ab_slice_fft = fftshift(fft2(Ex_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        Ey_ab_slice_fft = fftshift(fft2(Ey_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        Ez_ab_slice_fft = fftshift(fft2(Ez_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        
        E_0_film_slice_fft = fftshift(fft2(E_0_film_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        Ex_0_film_slice_fft = fftshift(fft2(Ex_0_film_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        Ey_0_film_slice_fft = fftshift(fft2(Ey_0_film_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        Ez_0_film_slice_fft = fftshift(fft2(Ez_0_film_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
        
        steady_z_diff = abs(Ez_ab(end-1, :, :)-temp_steady_z);
        steady_y_diff = abs(Ey_ab(end-1, :, :)-temp_steady_y);
        steady_z = [steady_z, sum(steady_z_diff(:))]; %#ok<AGROW>
        steady_y = [steady_y, sum(steady_y_diff(:))]; %#ok<AGROW>
        
        
        if save_dnymics_flag == 1
            display('saving field amplitude and N,Tc etc...');
            save([working_dictionary, '\output_data\E_0_film_', num2str(counter_I_gen), '.mat'], 'E_0_film');
            save([working_dictionary, '\output_data\Ex_0_film_', num2str(counter_I_gen), '.mat'], 'Ex_0_film');
            save([working_dictionary, '\output_data\Ey_0_film_', num2str(counter_I_gen), '.mat'], 'Ey_0_film');
            save([working_dictionary, '\output_data\Ez_0_film_', num2str(counter_I_gen), '.mat'], 'Ez_0_film');
            save([working_dictionary, '\output_data\F_gen_film.mat'], 'F_gen_film');
            save([working_dictionary, '\output_data\N_', num2str(counter_I_gen), '.mat'], 'N');
            save([working_dictionary, '\output_data\Tc_', num2str(counter_I_gen), '.mat'], 'Tc');
            save([working_dictionary, '\output_data\eps_', num2str(counter_I_gen), '.mat'], 'epsilon_excited_silicon');
        end
        
    end
    
    if see_fields_flag == 1 && mod(COUNTER_T, n_delta_t_av) == 0
        %if see_fields_flag==1
        %subplot(2,2,1);pcolor(reshape(double(Ex_SI(x_end_film-see_depth_animation,:,:)),N_points_y_all,N_points_z_all));colormap(jet(2^10));shading interp;axis equal;
        %subplot(2,2,2);pcolor(reshape(double(Ey_SI(x_end_film-see_depth_animation,:,:)),N_points_y_all,N_points_z_all));colormap(jet(2^10));shading interp;axis equal;
        %subplot(2,2,3);pcolor(reshape(double(Ez_SI(x_end_film-see_depth_animation,:,:)),N_points_y_all,N_points_z_all));colormap(jet(2^10));shading interp;axis equal;
        
        %subplot(2,2,4);pcolor(reshape(double(Ez_SI(x_start_film+(N_points_x_film+1-see_depth_animation):N_points_x_all,round(N_points_y_all/2),:)),...
        %    length(x_start_film+(N_points_x_film+1-see_depth_animation):N_points_x_all),N_points_z_all));colormap(jet(2^10));shading interp;
        %subplot(2,3,5);pcolor(reshape(double(Ex_SI(x_start_film+(N_points_x_film+1-see_depth_animation):N_points_x_all,round(N_points_y_all/2),:)),...
        %    length(x_start_film+(N_points_x_film+1-see_depth_animation):N_points_x_all),N_points_z_all));colormap(jet(2^10));shading interp;
        
        %subplot(2,3,1);pcolor(reshape(double(Ex_SI(x_end_film-see_depth_animation,N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge),...
        %    N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge))),...
        %    N_points_y_all-(2*(N_points_PML+1)+2*source_edge),N_points_z_all-(2*(N_points_PML+1)+2*source_edge)));colormap(jet(2^10));shading interp;axis equal;
        %subplot(2,3,2);pcolor(reshape(double(Ey_SI(x_end_film-see_depth_animation,N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge),...
        %    N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge))),...
        %    N_points_y_all-(2*(N_points_PML+1)+2*source_edge),N_points_z_all-(2*(N_points_PML+1)+2*source_edge)));colormap(jet(2^10));shading interp;axis equal;
        %subplot(2,3,3);pcolor(reshape(double(Ez_SI(x_end_film-see_depth_animation,N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge),...
        %    N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge))),...
        %    N_points_y_all-(2*(N_points_PML+1)+2*source_edge),N_points_z_all-(2*(N_points_PML+1)+2*source_edge)));colormap(jet(2^10));shading interp;axis equal;
        
        subaxis(2, 5, 1, 'spacing', 0, 'Margin', 0); pcolor(reshape(double(Ex_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge)), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge)); colormap(jet(2^10)); shading interp; axis off;%axis equal;
        subaxis(2, 5, 2, 'spacing', 0, 'Margin', 0); pcolor(reshape(double(Ey_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge)), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge)); colormap(jet(2^10)); shading interp; axis off;%axis equal;
        subaxis(2, 5, 3, 'spacing', 0, 'Margin', 0); pcolor(reshape(double(Ez_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge)), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge)); colormap(jet(2^10)); shading interp; axis off;%axis equal;
        
        subaxis(2, 5, 4, 'spacing', 0, 'Margin', 0); pcolor(reshape(double(E_0_film(end-see_depth_animation, source_edge+1:end-source_edge, source_edge+1:end-source_edge)), ...
            N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge)); colormap(jet(2^10)); shading interp; axis off;%axis equal;
        
        subaxis(2, 5, 5, 'spacing', 0, 'Margin', 0); pcolor(rough_height); shading interp; colormap(jet(2^10)); axis off;
        
        subaxis(2, 5, 6, 'spacing', 0, 'Margin', 0); pcolor(F_y, F_x, ((abs(double(Ex_0_film_slice_fft))))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)');
        caxis([0, color_factor * max(max(abs(double(Ex_0_film_slice_fft(colormap_F)))))]); axis off; xlim([-5, 5]); ylim([-5, 5]);
        %ylabel('k_y (k_0)');axis equal;
        %set(gca,'XTick',-7:1:7,'YTick',-7:1:7);
        
        subaxis(2, 5, 7, 'spacing', 0, 'Margin', 0); pcolor(F_y, F_x, ((abs(double(Ey_0_film_slice_fft))))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)');
        caxis([0, color_factor * max(max(abs(double(Ey_0_film_slice_fft(colormap_F)))))]); axis off; xlim([-5, 5]); ylim([-5, 5]);
        %ylabel('k_y (k_0)');axis equal;
        %set(gca,'XTick',-7:1:7,'YTick',-7:1:7);
        
        subaxis(2, 5, 8, 'spacing', 0, 'Margin', 0); pcolor(F_y, F_x, ((abs(double(Ez_0_film_slice_fft))))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)');
        caxis([0, color_factor * max(max(abs(double(Ez_0_film_slice_fft(colormap_F)))))]); axis off; xlim([-5, 5]); ylim([-5, 5]);
        %ylabel('k_y (k_0)');axis equal;
        %set(gca,'XTick',-7:1:7,'YTick',-7:1:7);
        
        subaxis(2, 5, 9, 'spacing', 0, 'Margin', 0); pcolor(F_y, F_x, ((abs(double(E_0_film_slice_fft))))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)');
        caxis([0, color_factor * max(max(abs(double(E_0_film_slice_fft(colormap_F)))))]); axis off; xlim([-5, 5]); ylim([-5, 5]);
        %ylabel('k_y (k_0)');axis equal;
        %set(gca,'XTick',-7:1:7,'YTick',-7:1:7);
        
        subaxis(2, 5, 10, 'spacing', 0, 'Margin', 0);
        pcolor(F_y, F_x, (abs(double(rough_height_fft)))); shading interp; colormap(jet(2^10));
        caxis([0, color_factor * max(max(abs(double(rough_height_fft(colormap_F)))))]); axis off; xlim([-5, 5]); ylim([-5, 5]);
        
        %getframe_nosteal_focus(h_animation);
        getframe(h_animation);
        if write_avi_flag == 1
            avi_obj = addframe(avi_obj, h_animation);
        end
        %getframe(h_animation);
        %caxis([min(min(Ez(x_end_film-1,:,:)))/10,max(max(Ez(x_end_film-1,:,:)))/10]);
        %figure(10);pcolor(reshape(double(Ey(:,:,round(N_points_z_all/2))),N_points_x_all,N_points_y_all));colormap(jet(2^10));shading interp;
        %caxis([min(min(Ey(:,:,round(N_points_z_all/2))))/10,max(max(Ey(:,:,round(N_points_z_all/2))))/10]);
        
        %subplot(1,3,1);pcolor(reshape(double(Ez_SI(:,:,round(N_points_z_all/2))),N_points_x_all,N_points_y_all));colormap(jet(2^10));shading interp;
        %subplot(1,3,2);pcolor(reshape(double(Ex_SI(:,:,round(N_points_z_all/2))),N_points_x_all,N_points_y_all));colormap(jet(2^10));shading interp;
        %subplot(1,3,3);pcolor(reshape(double(Ey_SI(:,:,round(N_points_z_all/2))),N_points_x_all,N_points_y_all));colormap(jet(2^10));shading interp;
        
        %caxis([min(min(gather(Ez(:,:,round(N_points_z_all/2)))))/10,max(max(gather(Ez(:,:,round(N_points_z_all/2)))))/10]);
        %caxis([-1e7 1e7]);
        %title('E_z');
        %caxis([-5e6 5e6]);
        %figure(10);pcolor(reshape(Dz(source_pos_x_center,:,:),N_points_y_all,N_points_z_all));caxis([-max(A(:)),max(A(:))]);shading interp;axis equal
        %getframe(h_animation);
    end
    
    if mod(COUNTER_T, round(1*delta_counter_t_recording)) == 0
        display(sprintf('%f %%', 100*(COUNTER_T / t_steps_total)));
        if record_field_flag == 1 && COUNTER_T >= t_steps_total_pulse - n_delta_t_av && COUNTER_T <= t_steps_total_pulse
            t_points_record = t_points_record + 1;
            Ez_SI_record_xyzt(:, :, :, t_points_record) = reshape((Ez_SI(N_points_PML+1+1:end-(N_points_PML + 1), N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge), ...
                N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge))), ...
                N_points_x_all-(2 * (N_points_PML + 1)), N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge));
            Ex_SI_record_xyzt(:, :, :, t_points_record) = reshape((Ex_SI(N_points_PML+1+1:end-(N_points_PML + 1), N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge), ...
                N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge))), ...
                N_points_x_all-(2 * (N_points_PML + 1)), N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge));
            Ey_SI_record_xyzt(:, :, :, t_points_record) = reshape((Ey_SI(N_points_PML+1+1:end-(N_points_PML + 1), N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge), ...
                N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge))), ...
                N_points_x_all-(2 * (N_points_PML + 1)), N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge));
            %Ez_SI_record_xyt(:,:,t_points_record)=reshape((Ez_SI(N_points_PML+1+1:end-(N_points_PML+1),N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge),round(N_points_z_all/2))),...
            %    N_points_x_all-(2*(N_points_PML+1)),N_points_y_all-(2*(N_points_PML+1)+2*source_edge));
            %Ez_SI_record_xzt(:,:,t_points_record)=reshape((Ez_SI(N_points_PML+1+1:end-(N_points_PML+1),round(N_points_y_all/2),N_points_PML+1+1+source_edge:end-(N_points_PML+1+source_edge))),...
            %    N_points_x_all-(2*(N_points_PML+1)),N_points_z_all-(2*(N_points_PML+1)+2*source_edge));
        end
    end
    
    t = t + delta_t; COUNTER_T = COUNTER_T + 1;
    FDTD_speed_inside = 1e-6 * (N_points_x_all * N_points_y_all * N_points_z_all) / time_used_inside;
    display(sprintf('%d, %f Mcells/s', COUNTER_T, FDTD_speed_inside));
    pause(0.1);
    
end

if write_avi_flag == 1
    avi_obj = close(avi_obj);
end
end


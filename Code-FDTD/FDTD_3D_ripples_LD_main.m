%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%3D-FDTD code (main file, put material properties, laser specifications etc.. in this code)
%%%to simulate laser pulse propagation in a lossy and
%%%dispersive medium using the Drude model, the transient permittivity change
%%%during the laser pulse is possible to be added in with a dynamically
%%%changing material properties.
%%%Author: Dr.Hao Zhang
%%%Date: March 17, 2014
%%% Lorentz-Drude model added, April 17, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all

%% Parallel pool checks
parallel_computing = true; % Variable encoding the will for parallel computing

if parallel_computing
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    if poolsize == 0
        parpool('local', 4);
    end
end

global COUNTER_T; %#ok<NUSED>

%material='W'; %%% selecting material
material = 'Ni-25000K';
order_Snm = 3; %%% the order of the taloy experssion used in the Snm function
display(['simulation material: ', material]);
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

working_dictionary = '.'
margins_ = ((1920 * 1080) / (1920 * 1080)) * [-2, -1.8, -0.5, 0.1]; papersize_ = ((1920 * 1080) / (1920 * 1080)) * [20, 20 / 2^0.5]; font_size_ = 20; marker_size_ = 6; line_width_ = 0.5;
see_fields_flag = 1; %%% see the Ex,Ey,Ez,Ex0,Ey0,Ez0 during the calculation or not
substracting_source_flag = 0; %%% substrating the source in the fft
write_avi_flag = 0; %%% wirte the animation to avi file or not
gpu_flag = 0; %%% 1-use gpu to do FDTD, 0-use CPU to do FDTD
Drude_LD_mode = 2; %%% 1-use Drude model (sigma,Chi,t0 form), 2-use Drude (the same form as lorentz-drude but with omega0=0) 3-use lorentz-drude model with 5 poles
single_precision_flag = 1; %%% 1-use single precision 0-use double precision
add_source = 1; %%% 1-add source 0-do not add source
refreshing_flag = 0; %%% 1-dielectric dynamics 0-no dielectric dynamics

pulse_number = 1; %%% applied fs laser pulse number

%multi_mode_type=1; %%% Eab2surface_morphology
multi_mode_type = 2; %%% Eab2rough
continue_mode = 0; %%% continue from the morpholoy of the crashed pulse
continue_pulse_num = 2; %%% continue from which pulse
ablation_depth = 1; %%% measured in delta_x, in order to get E_threshold, at E_ab=E_threshold, the ablation depth without surface roughness is the parameter ablation_depth
diffusion_flag = 0; %%% 1-with carrer diffusion, temperature diffusion etc... 0- without these diffusion
add_roughness_flag = 1; %%% 1-add surface roughness 0-no surface roughness
if strcmp(material, 'Vacuum') == 1
    add_roughness_flag = 0;
end
if refreshing_flag == 1
    Drude_LD_mode = 1;
end

roughness_type = 1; %%% 1-random roughness 2-gratings 3-Gaussian 4-Read image


%mean_rough=0.0001;%0.9999; %%% the mean of the roughness function
mean_rough = 0.2; %0.2; %%% the mean of the roughness function
rough_size = 1; %%% the size of a single roughness, measured in delta_x

%%% On considere 2 cellules pour créer la surface rugueuse (rough_thick).

rough_thick = 2; %%% the thickness of the roughness layer, including the grating height, measured in delta_x
base_thick = 8; %%% in the film where there is no roughness (for the first pulse)
prepare_to_ablate_thick = 2; %%% curcial for 1 pulse before the recyle use of the film
base_thick = base_thick + prepare_to_ablate_thick;
grating_roughness_thick = 1; %%% the thickness of the random roughness on the top of the grating, measured in delta_x, 0 is no random roughness
grating_period = 600e-9; %%% grating period
grating_width = 450e-9; %%% the width of the grating
grating_direction = 1; %%% 1- along y, 2- along z
load_roughness_flag = 0; %%% 1-use loaded roughness function 0-use random roughnes
roughness_dic = [working_dictionary, '/save_temp/']; %%% dictionary of the roughness function which is going to be loaded
save_dnymics_flag = 0; %%% 1-save E0, N, Tc etc... once half an optical cycle, 0-do not save
F_av = 0.2; %%% incident fluence J/cm^2
Tao = 10e-15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Wavelength %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 800e-9; %%% laser wavelength m %%%cc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

laser_calculation_time = 3 * Tao;
laser_calculation_time = 10e-15;
T_simulation = Tao / 2;
%T_simulation=2e-15;
%see_depth=0;
if pulse_number > 1 && multi_mode_type == 1
    rough_thick_multi_mode = 8;
else
    rough_thick_multi_mode = 0;
end
where_to_look = 0; %%% where to look at in the depth, +0 is immediately below the roughness, it is also the slice used to generate surface morphology used for the 2nd pulse in the multi-shots mode
see_depth = rough_thick + where_to_look;
%see_depth_animation=rough_thick+5;
see_depth_animation = see_depth;
%T_av=0.5;%%% 0.5--averge over half an optical cycle, 1- average over one optical cycle, it is also the length of the record_field
T_av = 1; %%% 0.5--averge over half an optical cycle, 1- average over one optical cycle, it is also the length of the record_field
record_field_flag = 1; %%% record the field of the last half optical cycle at a certain depth, very memory consuming.
%frequency_field_recording=5; %%% record the field every 1/frequency_field_recording optical cycle, if Inf, then every FDTD time step
frequency_field_recording = Inf; %%% record the field every 1/frequency_field_recording optical cycle, if Inf, then every FDTD time step
position_recording = see_depth_animation; %%% where to record in the depth, +0 is immediately below the roughness.
Besure_flag = 1; %%% 1--Besure the film thickness is correct. 2--Besure the SiO2 thickness is correct.
%delta_x_aiming=15e-9; %%% m the aimed delta_x
%delta_x_aiming=25e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Resolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_x_aiming = 40e-9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

air_thickness = 3 * delta_x_aiming; %%% m
film_thickness = (rough_thick + rough_thick_multi_mode + base_thick) * delta_x_aiming; %%% m
SiO2_thickness = 0 * delta_x_aiming; %%% m
Substract_thickness = 1 * delta_x_aiming; %%% m

%air_thickness=delta_x_aiming; %%% m
%film_thickness=delta_x_aiming; %%% m
%SiO2_thickness=delta_x_aiming; %%% m
%Substract_thickness=delta_x_aiming; %%% m

if Besure_flag == 1
    N_points_x_film_aiming = round(film_thickness/delta_x_aiming);
    delta_x = film_thickness / N_points_x_film_aiming; %%% the actual delta_x
elseif Besure_flag == 2
    N_points_x_SiO2_aiming = round(SiO2_thickness/delta_x_aiming);
    delta_x = SiO2_thickness / N_points_x_SiO2_aiming; %%% the actual delta_x
end

PML_thickness = 7 * delta_x; %% m
PML_coeff = 0.333; %%% xn=PML_coeff*((1:N_points_PML+1)/(N_points_PML+1)).^PML_order;
PML_order = 3; %%% xn=PML_coeff*((1:N_points_PML+1)/(N_points_PML+1)).^PML_order;

points_no_roughness_edge = 0e-9; %%% unit in m
no_rough_material = 0; %%% 1-the material of the no rough region is the simuated material, 0-the material of the no rough region is air
%sor_air_edge=200e-9; %%% srounding the structure with air before the PML layers
sor_air_edge = 0;
source_edge = max(sor_air_edge, points_no_roughness_edge) + 200e-9; %%% for the soft source, the distance between the edge of the source and the start of the PML layer, measured in delta_x

points_no_roughness_edge = round(points_no_roughness_edge/delta_x); %%% unit in delta_x
sor_air_edge = round(sor_air_edge/delta_x); %%% sorunding the structure with air before the PML layers
source_edge = round(source_edge/delta_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Roughness length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_dimension_main = 12e-6; %%% the size of the computation window in Y %12e-6
Z_dimension_main = 12e-6; %%% the size of the computation window in Z %12e-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y_dimension_main = Y_dimension_main + 2 * source_edge * delta_x; %%% the size of the computation window in Y
Z_dimension_main = Z_dimension_main + 2 * source_edge * delta_x; %%% the size of the computation window in Z
D_focal = (Y_dimension_main - 2 * source_edge * delta_x) * 1 / 3;
%D_focal=6e-6;
gau_ord = 7; %%% super gaussian distribution of order gau_ord, for gau_ord=2, it is a gaussian distribution.
%criterion_Eab2rough2=1; if 1, will not remove any layer
raduis_Eab2rough2 = 0.39; %%% if the raduis of continuous zeros is larger than this value times the total raduis, delete this layer and add another layer with all 1 to the bottom, used in the Eab2rough2 function, memory saving purpose.
w0 = D_focal / (4^(1 / gau_ord)); %%% beam waist at the focus

n_2_cm = 0.5e-14; %%% cm^2/W optical Kerr coefficient
n_2 = 1e-4 * n_2_cm; %%% m^2/W optical Kerr coefficient
beta_cm = 1.85e-9; %%% cm/W
beta = 1e-2 * beta_cm; %%% m/W
theta_impact_cm = 20;
theta_impact = 1e-4 * theta_impact_cm; %%% impact ionization coefficient m^2/J, the rate is theta*I
Tao_c_l_0 = 10e-12;
N_cri_phonon = 1e50; %%% m^-3
if Drude_LD_mode == 3 %%% Lorentz-Drude model parameters
    %epsilon_r_silicon=13.6+0.048i;%%% dielectric constant of c-silicon at 800nm
    epsilon_r_silicon = 1;
    epsilon_r_air = 1;
    %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
    %epsilon_r_glass=1;
    epsilon_r_glass = epsilon_r_silicon;
    omega0 = 0;
    if strcmp(material, 'Au') == 1
        %%% Lorentz-Drude parameters %%%
        %%% gold Au%%%
        %%% Drude
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        omegap = 9.03; %%% eV
        omegap = (omegap * eV) / hbar; %%% rad/s
        N_initial = (omegap^2 * m_ * m_e * epsilon_0) / couloum_e^2; %%%intrinsic carrier number m^-3
        gama0 = 0.053; %%% eV
        gama0 = (gama0 * eV) / hbar; %%%  rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
        tao_d = 1 / gama0; %%% collision time s/rad
        f0 = 0.760;
        %%% Lorentz
        f1 = 0.024;
        %f1=0;
        gama1 = 0.241; %%% eV
        omega1 = 0.415; %%% eV
        f2 = 0.01;
        %f2=0;
        gama2 = 0.345; %%% eV
        omega2 = 0.83; %%% eV
        f3 = 0.071;
        %f3=0;
        gama3 = 0.870; %%% eV
        omega3 = 2.969; %%% eV
        f4 = 0.601;
        %f4=0;
        gama4 = 2.494; %%% eV
        omega4 = 4.304; %%% eV
        f5 = 4.384;
        %f5=0;
        gama5 = 2.214; %%% eV
        omega5 = 13.32; %%% eV
    elseif strcmp(material, 'W')
        %%% Lorentz-Drude parameters %%%
        %%% Tungsten W %%%
        %%% Drude
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        omegap = 13.22; %%% eV
        omegap = (omegap * eV) / hbar; %%% rad/s
        N_initial = (omegap^2 * m_ * m_e * epsilon_0) / couloum_e^2; %%%intrinsic carrier number m^-3
        gama0 = 0.064; %%% eV
        gama0 = (gama0 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
        tao_d = 1 / gama0; %%% collision time s/rad
        f0 = 0.206;
        %%% Lorentz
        f1 = 0.054;
        gama1 = 0.530; %%% eV
        omega1 = 1.004; %%% eV
        f2 = 0.166;
        gama2 = 1.281; %%% eV
        omega2 = 1.917; %%% eV
        f3 = 0.706;
        gama3 = 3.332; %%% eV
        omega3 = 3.580; %%% eV
        f4 = 2.590;
        %f4=0;
        gama4 = 5.836; %%% eV
        omega4 = 7.498; %%% eV
        f5 = 0;
        gama5 = 0; %%% eV
        omega5 = 0; %%% eV
    elseif strcmp(material, 'W-25000K')
        %%% Lorentz-Drude parameters %%%
        %%% Tungsten W %%%
        %%% Drude
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        omegap = 12.6963; %%% eV
        omegap = (omegap * eV) / hbar; %%% rad/s
        N_initial = (omegap^2 * m_ * m_e * epsilon_0) / couloum_e^2; %%%intrinsic carrier number m^-3
        gama0 = 0.2440; %%% eV
        gama0 = (gama0 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
        tao_d = 1 / gama0; %%% collision time s/rad
        f0 = 0.3367;
        %%% Lorentz
        f1 = 0.2226;
        gama1 = 2.2411; %%% eV
        omega1 = 1.3349; %%% eV
        f2 = 0.0654;
        gama2 = 0.8453; %%% eV
        omega2 = 2.3722; %%% eV
        f3 = 0.2475;
        gama3 = 2.0119; %%% eV
        omega3 = 3.4374; %%% eV
        f4 = 1.2375;
        %f4=0;
        gama4 = 2.7033; %%% eV
        omega4 = 5.7975; %%% eV
        f5 = 0;
        gama5 = 0; %%% eV
        omega5 = 0; %%% eV
    elseif strcmp(material, 'Ni')
        %%% Lorentz-Drude parameters %%%
        %%% Nickel Ni%%%
        %%% Drude
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        omegap = 15.92; %%% eV
        %omegap=9;
        omegap = (omegap * eV) / hbar; %%% rad/s
        N_initial = (omegap^2 * m_ * m_e * epsilon_0) / couloum_e^2; %%%intrinsic carrier number m^-3
        gama0 = 0.048; %%% eV
        %gama0=0.053;
        gama0 = (gama0 * eV) / hbar; %%%  rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
        tao_d = 1 / gama0; %%% collision time s/rad
        %f0=0.76;
        f0 = 0.096;
        %%% Lorentz
        f1 = 0.100;
        %f1=0;
        gama1 = 4.511; %%% eV
        %gama1=1;
        omega1 = 0.174; %%% eV
        f2 = 0.135;
        %f2=0;
        gama2 = 1.334; %%% eV
        omega2 = 0.582; %%% eV
        f3 = 0.106;
        %f3=0;
        gama3 = 2.178; %%% eV
        omega3 = 1.597; %%% eV
        f4 = 0.729;
        %f4=0;
        gama4 = 6.292; %%% eV
        omega4 = 6.089; %%% eV
        f5 = 0;
        gama5 = 0; %%% eV
        omega5 = 0; %%% eV
    elseif strcmp(material, 'Ni-25000K') == 1
        %%% Lorentz-Drude parameters %%%
        %%% Nickel Ni-25000K%%%
        %%% Drude
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        omegap = 12.6711; %%% eV
        omegap = (omegap * eV) / hbar; %%% rad/s
        N_initial = (omegap^2 * m_ * m_e * epsilon_0) / couloum_e^2; %%%intrinsic carrier number m^-3
        gama0 = 0.4011; %%% eV
        gama0 = (gama0 * eV) / hbar; %%%  rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
        tao_d = 1 / gama0; %%% collision time s/rad
        f0 = 0.2344;
        %%% Lorentz
        f1 = 0.0089;
        gama1 = 0.1461; %%% eV
        omega1 = 1.4594; %%% eV
        f2 = 0.0061;
        gama2 = 0.1407; %%% eV
        omega2 = 1.9092; %%% eV
        f3 = 0.4122;
        gama3 = 3.7295; %%% eV
        omega3 = 2.5730; %%% eV
        f4 = 0.9909;
        gama4 = 6.0496; %%% eV
        omega4 = 6.3676; %%% eV
        f5 = 0;
        gama5 = 0; %%% eV
        omega5 = 0; %%% eV
    elseif strcmp(material, 'Vacuum') == 1
        %%% Lorentz-Drude parameters %%%
        %%% Nickel Ni-25000K%%%
        %%% Drude
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        omegap = 0; %%% eV
        omegap = (omegap * eV) / hbar; %%% rad/s
        N_initial = (omegap^2 * m_ * m_e * epsilon_0) / couloum_e^2; %%%intrinsic carrier number m^-3
        gama0 = 0.4011; %%% eV
        gama0 = (gama0 * eV) / hbar; %%%  rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
        tao_d = 1 / gama0; %%% collision time s/rad
        f0 = 1;
        %%% Lorentz
        f1 = 0;
        gama1 = 0; %%% eV
        omega1 = 0; %%% eV
        f2 = 0;
        gama2 = 0; %%% eV
        omega2 = 0; %%% eV
        f3 = 0;
        gama3 = 0; %%% eV
        omega3 = 0; %%% eV
        f4 = 0;
        gama4 = 0; %%% eV
        omega4 = 0; %%% eV
        f5 = 0;
        gama5 = 0; %%% eV
        omega5 = 0; %%% eV
    end
    gama1 = (gama1 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
    omega1 = (omega1 * eV) / hbar; %%% rad/s
    gama2 = (gama2 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
    omega2 = (omega2 * eV) / hbar; %%% rad/s
    gama3 = (gama3 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
    omega3 = (omega3 * eV) / hbar; %%% rad/s
    gama4 = (gama4 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
    omega4 = (omega4 * eV) / hbar; %%% rad/s
    gama5 = (gama5 * eV) / hbar; %%% rad/s %%% The paper (Aleksandar D. etal.) uses weird unit, they never use h, they use hbar for all cases
    omega5 = (omega5 * eV) / hbar; %%% rad/s
elseif Drude_LD_mode == 1 || 2 %%% Drude model parameters
    f0 = 1; %%% for Drude_LD_mode==2, f0 always equal 1
    omega0 = 0; %%% for Drude_LD_mode==2, omega0 always equal 0
    if strcmp(material, 'Si') == 1 %%% use it as Drude model
        %%% Drude
        epsilon_r_silicon = 13.6 + 0.048i; %%% dielectric constant of c-silicon at 800nm
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 0.18; %%% optical effective mass m_=m_1+m_k*Tc;
        %m_1=0.15;
        %m_k=3.1e-5;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 1e12; %%%intrinsic carrier number m^-3
        %N_initial=10e27;
        tao_d = 1.1e-15; %%% collision time s
    elseif strcmp(material, 'Au') == 1 %%% use it as Drude model
        %%% Drude
        epsilon_r_silicon = 7.6; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 5.9108e28; %%%intrinsic carrier number m^-3
        %N_initial=0;
        tao_d = 6.9334e-15; %%% collision time s
    elseif strcmp(material, 'W') == 1 %%% use it as Drude model
        %%% Drude
        epsilon_r_silicon = 11.35; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 1.2688e+29; %%%intrinsic carrier number m^-3
        tao_d = 1.2727e-16; %%% collision time s
    elseif strcmp(material, 'W_model') == 1 %%% use it as Drude model
        %%% Drude
        epsilon_r_silicon = 247.2327; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 0.3295; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 1.3983e+29; %%%intrinsic carrier number m^-3
        tao_d = 5.1280e-15; %%% collision time s
    elseif strcmp(material, 'W_25000K_model') == 1 %%% use it as Drude model
        %%% Drude
        epsilon_r_silicon = 145.0869; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 0.5516; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 1.5114e29; %%%intrinsic carrier number m^-3
        tao_d = 4.0323e-15; %%% collision time s
    elseif strcmp(material, 'Ni') == 1
        %%% Drude
        epsilon_r_silicon = 10; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 7.4107e28; %%%intrinsic carrier number m^-3
        tao_d = 4.6522e-16; %%% collision time s
    elseif strcmp(material, 'Ni-25000K') == 1
        %%% Drude
        epsilon_r_silicon = 10; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 5.0848e28; %%%intrinsic carrier number m^-3
        tao_d = 4.8708e-16; %%% collision time s
    elseif strcmp(material, '2+3i') == 1
        %%% Drude
        epsilon_r_silicon = 1; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 5.2268e28; %%%intrinsic carrier number m^-3
        tao_d = 2.1230e-16; %%% collision time s
    elseif strcmp(material, '3+2i') == 1
        %%% Drude
        epsilon_r_silicon = 10; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 5.8889e+28; %%%intrinsic carrier number m^-3
        tao_d = 1.7692e-16; %%% collision time s
    elseif strcmp(material, '3+0.1i') == 1
        %%% Drude
        epsilon_r_silicon = 10; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 2.3807e+27; %%%intrinsic carrier number m^-3
        tao_d = 7.1475e-16; %%% collision time s
    elseif strcmp(material, '3+1i') == 1
        %%% Drude
        epsilon_r_silicon = 10; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 3.4845e+28; %%%intrinsic carrier number m^-3
        tao_d = 1.4154e-16; %%% collision time s
    elseif strcmp(material, 'BMG') == 1
        %%% Drude
        epsilon_r_silicon = 10; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 6.2932e28; %%%intrinsic carrier number m^-3
        tao_d = 5.8497e-16; %%% collision time s
    elseif strcmp(material, 'Vacuum') == 1 %%% use it as Drude model
        %%% Drude
        epsilon_r_silicon = 1; %%% dielectric constant of c-silicon at 800nm, it is the epsilon infinity
        epsilon_r_air = 1;
        %epsilon_r_glass=1.453^2;%%%dielectric constant of fused silica at 800nm %%%
        %epsilon_r_glass=1;
        epsilon_r_glass = epsilon_r_silicon;
        m_1 = 1; %%% optical effective mass m_=m_1+m_k*Tc;
        m_k = 0;
        m_ = m_1 + m_k * 300;
        N_initial = 0; %%%intrinsic carrier number m^-3
        tao_d = 1e-15; %%% collision time s
    end
    f1 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    f2 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    f3 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    f4 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    f5 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    gama1 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    omega1 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    gama2 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    omega2 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    gama3 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    omega3 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    gama4 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    omega4 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    gama5 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
    omega5 = 0; %%% for Drude model, those parameters are never used in the loop function, but still has a value to ensure the loop function reveive enough input
end

save([working_dictionary, '/save_temp/delta_x.mat'], 'delta_x');
save([working_dictionary, '/save_temp/rough_thick.mat'], 'rough_thick');
save([working_dictionary, '/save_temp/points_no_roughness_edge.mat'], 'points_no_roughness_edge');
save([working_dictionary, '/save_temp/sor_air_edge.mat'], 'sor_air_edge');
save([working_dictionary, '/save_temp/source_edge.mat'], 'source_edge');
save([working_dictionary, '/save_temp/rough_thick_multi_mode.mat'], 'rough_thick_multi_mode');
f = c_vacuum / lambda;
lambda_range = linspace(0.8*lambda, 1.2*lambda, 1e3);
f_range = c_vacuum ./ lambda_range;
omegap = ((N_initial * couloum_e^2) ./ (m_ * m_e * epsilon_0)).^0.5;
if Drude_LD_mode == 3
    eps_exci = epsilon_r_silicon + (f0 * omegap.^2) ./ (omega0^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama0) + ...
        (f1 * omegap.^2) ./ (omega1^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama1) + ...
        (f2 * omegap.^2) ./ (omega2^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama2) + ...
        (f3 * omegap.^2) ./ (omega3^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama3) + ...
        (f4 * omegap.^2) ./ (omega4^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama4) + ...
        (f5 * omegap.^2) ./ (omega5^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama5);
    eps_exci_range = epsilon_r_silicon + (f0 * omegap.^2) ./ (omega0^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama0) + ...
        (f1 * omegap.^2) ./ (omega1^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama1) + ...
        (f2 * omegap.^2) ./ (omega2^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama2) + ...
        (f3 * omegap.^2) ./ (omega3^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama3) + ...
        (f4 * omegap.^2) ./ (omega4^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama4) + ...
        (f5 * omegap.^2) ./ (omega5^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama5);
elseif Drude_LD_mode == 1
    eps_exci = epsilon_r_silicon - (((omegap) ./ (2 * pi * f)).^2) .* (1 ./ (1 + 1i * (1 ./ (2 * pi * f * tao_d))));
    eps_exci_range = epsilon_r_silicon - (((omegap) ./ (2 * pi * f_range)).^2) .* (1 ./ (1 + 1i * (1 ./ (2 * pi * f_range * tao_d))));
elseif Drude_LD_mode == 2
    gama0 = 1 / tao_d;
    eps_exci = epsilon_r_silicon + (f0 * omegap.^2) ./ (omega0^2 - (2 * pi * f)^2 - 1i * 2 * pi * f * gama0); %%% delectric constant of excited silicon
    eps_exci_range = epsilon_r_silicon + (f0 * omegap.^2) ./ (omega0^2 - (2 * pi * f_range).^2 - 1i * 2 * pi * f_range * gama0); %%% delectric constant of excited silicon
end
display(sprintf('rough_size=%f nm', 1e9*rough_size*delta_x));
display(sprintf('rough_thick=%f nm', 1e9*rough_thick*delta_x));
display(sprintf('epsion_exci=%f+%fi', real(eps_exci), imag(eps_exci)));
figure('name', 'eps(omega)');
plot(lambda_range*1e9, real(eps_exci_range));
hold on;
plot(lambda_range*1e9, imag(eps_exci_range), '--');
xlabel('\lambda (nm)'); ylabel('\epsilon (\lambda)');

t = cputime;
if continue_mode == 1
    pulse_num_range = continue_pulse_num:pulse_number;
    load_roughness_flag = 1;
    roughness_function = importdata([working_dictionary, '/save_temp/roughness_function_pulse_', num2str(continue_pulse_num-1), '.mat']);
    save([working_dictionary, '/save_temp/roughness_function.mat'], 'roughness_function');
    E_threshold = importdata([working_dictionary, '/save_temp/E_threshold.mat']);
    N_points_x_film = importdata([working_dictionary, '/save_temp/N_points_x_film.mat']);
    for i = N_points_x_film + 1:-1:1
        if min(min(roughness_function(i, :, :))) == 1
            where_to_look = N_points_x_film + 1 - i - rough_thick - rough_thick_multi_mode; %%% search for the layer that is immediately below the new roughness
            break;
        end
    end
    see_depth = rough_thick + rough_thick_multi_mode + where_to_look;
    see_depth_animation = see_depth;
    position_recording = see_depth_animation;
else
    pulse_num_range = 0:pulse_number;
end

for pulse_num = pulse_num_range
    display(sprintf('The %d pulse', pulse_num));
    display(sprintf('continue_mode=%d', continue_mode));
    if pulse_num == 0
        add_roughness_flag = 0; %%% get the E_threshold using the simulation without roughness
    else
        add_roughness_flag = 1;
    end
    [Ez_SI_record_xyzt, Ex_SI_record_xyzt, Ey_SI_record_xyzt, steady_z, steady_y, E_ab, Ex_ab, Ey_ab, Ez_ab, E_0_film, Ex_0_film, Ey_0_film, Ez_0_film, Ex_SI, Ey_SI, Ez_SI, Hx, Hy, Hz, N, epsilon_r_all, air_thickness_act, film_thickness_act, SiO2_thickness_act, ...
        Substract_thickness_act, PML_thickness_act, N_points_x_air, roughness_function, x_corr, y_corr, z_corr, N_points_x_film, N_points_x_SiO2, ...
        N_points_x_Substract, N_points_PML, x_end_film, N_points_x_all, N_points_y_all, N_points_z_all, N_points_y_main, N_points_z_main, E_ab_t, E_Uc_t, E_Ul_t] = FDTD_3D_ripples_LD_loop ...
        (N_initial, m_, epsilon_r_silicon, epsilon_r_air, epsilon_r_glass, tao_d, f0, omega0, f1, gama1, omega1, f2, gama2, omega2, f3, gama3, omega3, f4, gama4, omega4, ...
        f5, gama5, omega5, order_Snm, air_thickness, film_thickness, SiO2_thickness, ...
        Substract_thickness, PML_thickness, PML_coeff, PML_order, Y_dimension_main, Z_dimension_main, delta_x, laser_calculation_time, T_simulation, add_source, source_edge, points_no_roughness_edge, no_rough_material, ...
        F_av, w0, gau_ord, Tao, lambda, add_roughness_flag, roughness_type, grating_roughness_thick, mean_rough, rough_size, rough_thick, grating_period, grating_width, grating_direction, load_roughness_flag, roughness_dic, ...
        m_1, m_k, Tao_c_l_0, N_cri_phonon, beta, theta_impact, n_2, refreshing_flag, gpu_flag, working_dictionary, ...
        save_dnymics_flag, single_precision_flag, diffusion_flag, pulse_num, Drude_LD_mode, see_depth_animation, see_fields_flag, write_avi_flag, sor_air_edge, record_field_flag, T_av, frequency_field_recording, position_recording);
    if pulse_num == 0
        E_threshold = 1e-6 * max(max(E_ab(end-ablation_depth, :, :))) %%% J/cm^3, get the E_threshold using the simulation without roughness
        save([working_dictionary, '/save_temp/E_threshold.mat'], 'E_threshold');
        %save([working_dictionary,'roughness_function_1stpulse'],'roughness_function');
    end
    display('generating surface morphology...');
    
    %roughness_function=Eab2rough(E_threshold,E_ab,Ez_ab,roughness_function);
    if min(min(roughness_function(prepare_to_ablate_thick, :, :))) == 0 %%% if the pulse has drilled a hole in the bottom layer
        roughness_function = Eab2rough2(E_threshold, E_ab, Ez_ab, roughness_function, raduis_Eab2rough2);
    else
        roughness_function = Eab2rough(E_threshold, E_ab, Ez_ab, roughness_function);
    end
    %roughness_function=Eab2surface_morphology(E_ab,see_depth,(rough_thick+rough_thick_multi_mode)*delta_x,base_thick*delta_x,N_points_x_film,N_points_y_main,N_points_z_main,delta_x,pulse_num);
    if pulse_num >= 2
        for i = N_points_x_film + 1:-1:1
            if min(min(roughness_function(i, :, :))) == 1
                where_to_look = N_points_x_film + 1 - i - rough_thick - rough_thick_multi_mode; %%% search for the layer that is immediately below the new roughness
                break;
            end
        end
    end
    see_depth = rough_thick + rough_thick_multi_mode + where_to_look;
    see_depth_animation = see_depth;
    position_recording = see_depth_animation;
    save([working_dictionary, '/save_temp/roughness_function.mat'], 'roughness_function');
    if pulse_num >= 1
        load_roughness_flag = 1;
    end
    save([working_dictionary, '/save_temp/E_ab_pulse_', num2str(pulse_num), '.mat'], 'E_ab');
    save([working_dictionary, '/save_temp/Ex_ab_pulse_', num2str(pulse_num), '.mat'], 'Ex_ab');
    save([working_dictionary, '/save_temp/Ey_ab_pulse_', num2str(pulse_num), '.mat'], 'Ey_ab');
    save([working_dictionary, '/save_temp/Ez_ab_pulse_', num2str(pulse_num), '.mat'], 'Ez_ab');
    save([working_dictionary, '/save_temp/E_0_film_pulse_', num2str(pulse_num), '.mat'], 'E_0_film');
    save([working_dictionary, '/save_temp/Ex_0_film_pulse_', num2str(pulse_num), '.mat'], 'Ex_0_film');
    save([working_dictionary, '/save_temp/Ey_0_film_pulse_', num2str(pulse_num), '.mat'], 'Ey_0_film');
    save([working_dictionary, '/save_temp/Ez_0_film_pulse_', num2str(pulse_num), '.mat'], 'Ez_0_film');
    save([working_dictionary, '/save_temp/roughness_function_pulse_', num2str(pulse_num), '.mat'], 'roughness_function');
    if record_field_flag == 1
        save([working_dictionary, '/save_temp/Ez_SI_record_xyzt_', num2str(pulse_num), '.mat'], 'Ez_SI_record_xyzt');
        save([working_dictionary, '/save_temp/Ex_SI_record_xyzt_', num2str(pulse_num), '.mat'], 'Ex_SI_record_xyzt');
        save([working_dictionary, '/save_temp/Ey_SI_record_xyzt_', num2str(pulse_num), '.mat'], 'Ey_SI_record_xyzt');
    end
end
cputime - t

%save([working_dictionary,'Ex_SI.mat'], 'Ex');
%save([working_dictionary,'Ey_SI.mat'], 'Ey');
%save([working_dictionary,'Ez_SI.mat'], 'Ez');
%save([working_dictionary,'Hx_SI.mat'], 'Hx');
%save([working_dictionary,'Hy_SI.mat'], 'Hy');
%save([working_dictionary,'Hz_SI.mat'], 'Hz');

save([working_dictionary, '/save_temp/E_ab.mat'], 'E_ab');
save([working_dictionary, '/save_temp/Ex_ab.mat'], 'Ex_ab');
save([working_dictionary, '/save_temp/Ey_ab.mat'], 'Ey_ab');
save([working_dictionary, '/save_temp/Ez_ab.mat'], 'Ez_ab');
save([working_dictionary, '/save_temp/Ex_0_film.mat'], 'Ex_0_film');
save([working_dictionary, '/save_temp/Ey_0_film.mat'], 'Ey_0_film');
save([working_dictionary, '/save_temp/Ez_0_film.mat'], 'Ez_0_film');

%%% fft %%%
Ez_ab = 1e-6 * double((Ez_ab)); %%% in J/cm^3
Ez_ab_slice = reshape(Ez_ab(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
E_ab = 1e-6 * double((E_ab)); %%% in J/cm^3
E_ab_slice = reshape(E_ab(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
Ex_ab = 1e-6 * double((Ex_ab)); %%% in J/cm^3
Ex_ab_slice = reshape(Ex_ab(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
Ey_ab = 1e-6 * double((Ey_ab)); %%% in J/cm^3
Ey_ab_slice = reshape(Ey_ab(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
Ex_0_film_squared = (double((Ex_0_film))).^2; %%% in V^2/m^2, the field envelope squared in the last optical half cycle
Ex_0_film_squared_slice = reshape(Ex_0_film_squared(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
Ey_0_film_squared = (double((Ey_0_film))).^2; %%% in V^2/m^2, the field envelope squared in the last optical half cycle
Ey_0_film_squared_slice = reshape(Ey_0_film_squared(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
Ez_0_film_squared = (double((Ez_0_film))).^2; %%% in V^2/m^2, the field envelope squared in the last optical half cycle
Ez_0_film_squared_slice = reshape(Ez_0_film_squared(end-see_depth, source_edge+1:end-source_edge, source_edge+1:end-source_edge), ...
    N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);

Ex_SI_slice = double(reshape((Ex_SI(x_end_film-see_depth_animation, N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge), ...
    N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge))), ...
    N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge)));%%% V/m
Ey_SI_slice = double(reshape((Ey_SI(x_end_film-see_depth_animation, N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge), ...
    N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge))), ...
    N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge)));%%% V/m
Ez_SI_slice = double(reshape((Ez_SI(x_end_film-see_depth_animation, N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge), ...
    N_points_PML+1+1+source_edge:end-(N_points_PML + 1 + source_edge))), ...
    N_points_y_all-(2 * (N_points_PML + 1) + 2 * source_edge), N_points_z_all-(2 * (N_points_PML + 1) + 2 * source_edge)));%%% V/m

if substracting_source_flag == 1
    %%% get the fourier components of the soruce
    YZ = meshgrid(y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)));
    amp_guess_E_ab = max(E_ab_slice(:));
    amp_guess_Ez_ab = max(Ez_ab_slice(:));
    y_shift_guess = y_corr(round(length(y_corr)/2));
    z_shift_guess = z_corr(round(length(z_corr)/2));
    error_fun_E_ab = @(x0) x0(1) * exp(-(2 * (abs((YZ - x0(2))).^gau_ord + abs((YZ' - x0(3))).^gau_ord))./(x0(4)^gau_ord)) - E_ab_slice;
    error_fun_Ez_ab = @(x0) x0(1) * exp(-(2 * (abs((YZ - x0(2))).^gau_ord + abs((YZ' - x0(3))).^gau_ord))./(x0(4)^gau_ord)) - Ez_ab_slice;
    %[coeff_E_ab,~,~,~]=lsqnonlin(error_fun_E_ab,[amp_guess_E_ab y_shift_guess z_shift_guess w0]);
    [coeff_Ez_ab, resnorm, residual, exitflag] = lsqnonlin(error_fun_Ez_ab, [amp_guess_Ez_ab, y_shift_guess, z_shift_guess, w0]);
    E_ab_slice_fitted = coeff_E_ab(1) * exp(-(2 * (abs((YZ - coeff_E_ab(2))).^gau_ord + abs((YZ' - coeff_E_ab(3))).^gau_ord))./(coeff_E_ab(4)^gau_ord));
    Ez_ab_slice_fitted = coeff_Ez_ab(1) * exp(-(2 * (abs((YZ - coeff_Ez_ab(2))).^gau_ord + abs((YZ' - coeff_Ez_ab(3))).^gau_ord))./(coeff_Ez_ab(4)^gau_ord));
    %%% get the fourier components of the soruce
else
    E_ab_slice_fitted = zeros(N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
    Ez_ab_slice_fitted = zeros(N_points_y_main+1-2*source_edge, N_points_z_main+1-2*source_edge);
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
color_factor = 0.5;
%%% use to defind the colormap for the fft image

%Ez_ab_slice_fft=fftshift(fft2(Ez_ab_slice,N_fft_x,N_fft_y)/((N_points_y_main+1)*(N_points_z_main+1))); %%% fouirer transform of the E_ab_slice
Ez_ab_slice_fft = fftshift(fft2(Ez_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the Ez_ab_slice
Ey_ab_slice_fft = fftshift(fft2(Ey_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the Ez_ab_slice
Ex_ab_slice_fft = fftshift(fft2(Ex_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the Ez_ab_slice
E_ab_slice_fft = fftshift(fft2(E_ab_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice
Ex_0_film_slice_fft = fftshift(fft2(Ex_0_film_squared_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the Ez_0_film_slice
Ey_0_film_slice_fft = fftshift(fft2(Ey_0_film_squared_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the Ez_0_film_slice
Ez_0_film_slice_fft = fftshift(fft2(Ez_0_film_squared_slice, N_fft_x, N_fft_y)); %%% fouirer transform of the Ez_0_film_slice
Ex_SI_slice_fft = fftshift(fft2(Ex_SI_slice, N_fft_x, N_fft_y));
Ey_SI_slice_fft = fftshift(fft2(Ey_SI_slice, N_fft_x, N_fft_y));
Ez_SI_slice_fft = fftshift(fft2(Ez_SI_slice, N_fft_x, N_fft_y));

Ez_ab_slice_fitted_fft = fftshift(fft2(Ez_ab_slice_fitted, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice_fitted
Ez_ab_slice_fft = Ez_ab_slice_fft - Ez_ab_slice_fitted_fft; %%% substracting the fourier components of the source
E_ab_slice_fitted_fft = fftshift(fft2(E_ab_slice_fitted, N_fft_x, N_fft_y)); %%% fouirer transform of the E_ab_slice_fitted
E_ab_slice_fft = E_ab_slice_fft - E_ab_slice_fitted_fft; %%% substracting the fourier components of the source

h_Ex_yz = figure('name', 'Ex_yz_SI');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), 1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), Ex_SI_slice);
shading interp;colormap(jet(2^10));axis equal;
xlabel('y (\mum)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('z (\mum)', 'Interpreter', 'tex', 'FontSize', font_size_); set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ey_yz = figure('name', 'Ey_yz_SI');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), 1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), Ey_SI_slice);
shading interp;colormap(jet(2^10));axis equal;
xlabel('y (\mum)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('z (\mum)', 'Interpreter', 'tex', 'FontSize', font_size_); set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ez_yz = figure('name', 'Ez_yz_SI');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), 1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), Ez_SI_slice);
shading interp;colormap(jet(2^10));axis equal;
xlabel('y (\mum)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('z (\mum)', 'Interpreter', 'tex', 'FontSize', font_size_); set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ex_fft = figure('name', 'fft_Ex');
pcolor(F_y, F_x, (abs(Ex_SI_slice_fft))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ex_SI_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_Ey_fft = figure('name', 'fft_Ey');
pcolor(F_y, F_x, (abs(Ey_SI_slice_fft))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ey_SI_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_Ez_fft = figure('name', 'fft_Ez');
pcolor(F_y, F_x, (abs(Ez_SI_slice_fft))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ez_SI_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

%h_Hz_yz=figure('name','Hz_yz_SI');
%pcolor(1e6*z_corr,1e6*y_corr,reshape(double(Hz(x_end_film,:,:)),N_points_y_all,N_points_z_all));shading interp;colormap(hot(2^10));axis equal;
%xlabel('y (\mum)','Interpreter','tex','FontSize',font_size_);ylabel('z (\mum)','Interpreter','tex','FontSize',font_size_);set(gca,'fontsize',font_size_,'FontAngle','italic');

%h_Hy_yz=figure('name','Hy_yz_SI');
%pcolor(1e6*z_corr,1e6*y_corr,reshape(double(Hy(x_end_film,:,:)),N_points_y_all,N_points_z_all));shading interp;colormap(hot(2^10));axis equal;
%xlabel('y (\mum)','Interpreter','tex','FontSize',font_size_);ylabel('z (\mum)','Interpreter','tex','FontSize',font_size_);set(gca,'fontsize',font_size_,'FontAngle','italic');

%h_Hx_yz=figure('name','Hx_yz_SI');
%pcolor(1e6*z_corr,1e6*y_corr,reshape(double(Hx(x_end_film,:,:)),N_points_y_all,N_points_z_all));shading interp;colormap(hot(2^10));axis equal;
%xlabel('y (\mum)','Interpreter','tex','FontSize',font_size_);ylabel('z (\mum)','Interpreter','tex','FontSize',font_size_);set(gca,'fontsize',font_size_,'FontAngle','italic');

h_E_ab = figure('name', 'E_ab');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(E_ab_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ez_ab = figure('name', 'Ez_ab');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(Ez_ab_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ez_ab_fft = figure('name', 'fft_Ez_ab');
pcolor(F_y, F_x, (abs(double(Ez_ab_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ez_ab_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_Ey_ab = figure('name', 'Ey_ab');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(Ey_ab_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ey_ab_fft = figure('name', 'fft_Ey_ab');
pcolor(F_y, F_x, (abs(double(Ey_ab_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ey_ab_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_Ex_ab = figure('name', 'Ex_ab');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(Ex_ab_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_Ex_ab_fft = figure('name', 'fft_Ex_ab');
pcolor(F_y, F_x, (abs(double(Ex_ab_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ex_ab_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_E_ab_fft = figure('name', 'fft_E_ab');
pcolor(F_y, F_x, (abs(double(E_ab_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(E_ab_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_E0_x = figure('name', 'E0x_squared');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(Ex_0_film_squared_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_E0_y = figure('name', 'E0y_squared');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(Ey_0_film_squared_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_E0_z = figure('name', 'E0z_squared');
pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML + 1 + source_edge)), ...
    1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML + 1 + source_edge)), ...
    double(Ez_0_film_squared_slice));
shading interp;colormap(jet(2^10));axis equal;
xlabel('z (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_); ylabel('y (\mu m)', 'Interpreter', 'tex', 'FontSize', font_size_);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic');

h_E0_x_fft = figure('name', 'fft_E0x_squared');
pcolor(F_y, F_x, (abs(double(Ex_0_film_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ex_0_film_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_E0_y_fft = figure('name', 'fft_E0y_squared');
pcolor(F_y, F_x, (abs(double(Ey_0_film_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ey_0_film_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

h_E0_z_fft = figure('name', 'fft_E0z_squared');
pcolor(F_y, F_x, (abs(double(Ez_0_film_slice_fft)))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
caxis([0, color_factor * max(max(abs(double(Ez_0_film_slice_fft(colormap_F)))))]);
ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); axis equal; xlim([-5, 5]); ylim([-5, 5]);
set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -5:1:5, 'YTick', -5:1:5);

%fillPage(h_Ez_yz, 'margins',margins_, 'papersize', papersize_);
%saveas(h_Ez_yz,[working_dictionary,'Ez_yz'],'fig');
%print(h_Ez_yz,'-dpdf','-r300',[working_dictionary,'Ez_yz']);

%fillPage(h_Ey_yz, 'margins',margins_, 'papersize', papersize_);
%saveas(h_Ey_yz,[working_dictionary,'Ey_yz'],'fig');
%print(h_Ey_yz,'-dpdf','-r300',[working_dictionary,'Ey_yz']);

%fillPage(h_Ex_yz, 'margins',margins_, 'papersize', papersize_);
%saveas(h_Ex_yz,[working_dictionary,'Ex_yz'],'fig');
%print(h_Ex_yz,'-dpdf','-r300',[working_dictionary,'Ex_yz']);

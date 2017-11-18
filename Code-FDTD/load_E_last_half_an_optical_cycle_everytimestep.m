clear;clc;close all;
%working_dictionary='/Users/ludovic/Documents/MATLAB/Projet Recherche et Innovation/Version_Efield_time';
working_dictionary='C:/Data/Theses_Postdocs/Stage_ludovic/Projet_RI_Chopineau/Code-FDTD'
exci='/save_temp/';
margins_=[-2 -3 -0.5 0.1];papersize_=[18 18/2^0.5];font_size_=20;marker_size_=6;line_width_=0.5;
load_pulse_num=1; %%% which pulse to load
E_xyz_flag=3; %%% 1-see Ex field, 2-see Ey field, 3-see Ez field;

%Ex_SI_record_yzt_vaccum=importdata([working_dictionary,exci,'Ex_SI_record_yzt_vaccum.mat']); %%% Ex field every time step of the last half an optical cycle
%Ey_SI_record_yzt_vaccum=importdata([working_dictionary,exci,'Ey_SI_record_yzt_vaccum.mat']); %%% Ey field every time step of the last half an optical cycle
%Ez_SI_record_yzt_vaccum=importdata([working_dictionary,exci,'Ez_SI_record_yzt_vaccum.mat']); %%% Ez field every time step of the last half an optical cycle

if E_xyz_flag==1
    Ex_SI_record_xyzt=importdata([working_dictionary,exci,'Ex_SI_record_xyzt_',num2str(load_pulse_num),'.mat']); %%% Ex field every time step of the last half an optical cycle
elseif E_xyz_flag==2
    Ey_SI_record_xyzt=importdata([working_dictionary,exci,'Ey_SI_record_xyzt_',num2str(load_pulse_num),'.mat']); %%% Ey field every time step of the last half an optical cycle
elseif E_xyz_flag==3
    Ez_SI_record_xyzt=importdata([working_dictionary,exci,'Ez_SI_record_xyzt_',num2str(load_pulse_num),'.mat']); %%% Ez field every time step of the last half an optical cycle
end

%Ex_SI_record_yzt_sca=Ex_SI_record_yzt-Ex_SI_record_yzt_vaccum; %%% Ex scattered field every time step of the last half an optical cycle
%Ey_SI_record_yzt_sca=Ey_SI_record_yzt-Ey_SI_record_yzt_vaccum; %%% Ey scattered field every time step of the last half an optical cycle
%Ez_SI_record_yzt_sca=Ez_SI_record_yzt-Ez_SI_record_yzt_vaccum; %%% Ez scattered field every time step of the last half an optical cycle

y_corr=importdata([working_dictionary,exci,'y_corr.mat']);
z_corr=importdata([working_dictionary,exci,'z_corr.mat']);
N_points_PML=importdata([working_dictionary,exci,'N_points_PML.mat']);
N_points_z_all=importdata([working_dictionary,exci,'N_points_z_all.mat']);
N_points_y_all=importdata([working_dictionary,exci,'N_points_y_all.mat']);
N_points_x_all=importdata([working_dictionary,exci,'N_points_x_all.mat']);
N_points_x_air=importdata([working_dictionary,exci,'N_points_x_air.mat']);
N_points_x_SiO2=importdata([working_dictionary,exci,'N_points_x_SiO2.mat']);
N_points_x_Substract=importdata([working_dictionary,exci,'N_points_x_Substract.mat']);
N_points_x_film=importdata([working_dictionary,exci,'N_points_x_film.mat']);
N_points_y_main=importdata([working_dictionary,exci,'N_points_y_main.mat']);
N_points_z_main=importdata([working_dictionary,exci,'N_points_z_main.mat']);
delta_x=importdata([working_dictionary,exci,'delta_x.mat']); %%% in meter
rough_thick=importdata([working_dictionary,exci,'rough_thick.mat']);

points_no_roughness_edge=importdata([working_dictionary,exci,'points_no_roughness_edge.mat']);
sor_air_edge=importdata([working_dictionary,exci,'sor_air_edge.mat']);
source_edge=importdata([working_dictionary,exci,'source_edge.mat']);
roughness_function=importdata([working_dictionary,exci,'roughness_function','_pulse_',num2str(load_pulse_num),'.mat']); %%% roughness function

x_start_Substract=1;
x_end_Substract=x_start_Substract+N_points_x_Substract;
x_start_SiO2=x_end_Substract+1;
x_end_SiO2=x_start_SiO2+N_points_x_SiO2;
x_start_film=x_end_SiO2+1;
x_end_film=x_start_film+N_points_x_film;
x_start_air=x_end_film+1;
x_end_air=x_start_air+N_points_x_air;

see_position1=x_start_air;
see_position2=x_end_film-1;
see_position3=x_end_film-10;
size_data=size(Ez_SI_record_xyzt);
timesteps=size_data(4);

if E_xyz_flag==1
    Ex_SI_record_position1_yzt=reshape(Ex_SI_record_xyzt(see_position1,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
    Ex_SI_record_position2_yzt=reshape(Ex_SI_record_xyzt(see_position2,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
    Ex_SI_record_position3_yzt=reshape(Ex_SI_record_xyzt(see_position3,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
elseif E_xyz_flag==2
    Ey_SI_record_position1_yzt=reshape(Ey_SI_record_xyzt(see_position1,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
    Ey_SI_record_position2_yzt=reshape(Ey_SI_record_xyzt(see_position2,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
    Ey_SI_record_position3_yzt=reshape(Ey_SI_record_xyzt(see_position3,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
elseif E_xyz_flag==3
    Ez_SI_record_position1_yzt=reshape(Ez_SI_record_xyzt(see_position1,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
    Ez_SI_record_position2_yzt=reshape(Ez_SI_record_xyzt(see_position2,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
    Ez_SI_record_position3_yzt=reshape(Ez_SI_record_xyzt(see_position3,:,:,:),N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge,timesteps);
end
%E_SI_record_position1_yzt=sqrt(Ex_SI_record_position1_yzt.^2+Ey_SI_record_position1_yzt.^2+Ez_SI_record_position1_yzt.^2);
%E_SI_record_position2_yzt=sqrt(Ex_SI_record_position2_yzt.^2+Ey_SI_record_position2_yzt.^2+Ez_SI_record_position2_yzt.^2);
%E_SI_record_position3_yzt=sqrt(Ex_SI_record_position3_yzt.^2+Ey_SI_record_position3_yzt.^2+Ez_SI_record_position3_yzt.^2);

%theta_record_position1_yzt=acosd(Ez_SI_record_position1_yzt./E_SI_record_position1_yzt); %%% the angle between the E field and the z axis
%theta_record_position2_yzt=acosd(Ez_SI_record_position2_yzt./E_SI_record_position2_yzt); %%% the angle between the E field and the z axis
%theta_record_position3_yzt=acosd(Ez_SI_record_position3_yzt./E_SI_record_position3_yzt); %%% the angle between the E field and the z axis

%%% premeter for fft and ifft %%%
N_fft_x=max(2^9,N_points_y_main+1-2*source_edge);
N_fft_y=max(2^9,N_points_z_main+1-2*source_edge);

lambda=800e-9;
k_0=2*pi/lambda;
F_spatial=2*pi*(1/(delta_x))/k_0; %%% sampling frequency (spatial) normliazed to k_0
F_x=linspace(-F_spatial/2,F_spatial/2,N_fft_x); %%% coordinate in k space %%%
F_y=linspace(-F_spatial/2,F_spatial/2,N_fft_y); %%% coordinate in k space %%%
F_xy_mesh=meshgrid(F_x,F_y);
k_range=5; %%% k0, the range of wave number to see

%%% use to defind the colormap for the fft image
F_colormap=1; %%% k_0
colormap_F=find(abs(F_xy_mesh)>F_colormap);
color_factor=0.5;
%%% use to defind the colormap for the fft image

Ez_SI_record_position1_yzt_fft=zeros(N_fft_x,N_fft_y,size_data(3));
Ez_SI_record_position2_yzt_fft=zeros(N_fft_x,N_fft_y,size_data(3));
Ez_SI_record_position3_yzt_fft=zeros(N_fft_x,N_fft_y,size_data(3));

for i=1:timesteps
    Ez_SI_record_position1_yzt_fft(:,:,i)=fftshift(fft2(Ez_SI_record_position1_yzt(:,:,i),N_fft_x,N_fft_y));
    Ez_SI_record_position2_yzt_fft(:,:,i)=fftshift(fft2(Ez_SI_record_position2_yzt(:,:,i),N_fft_x,N_fft_y));
    Ez_SI_record_position3_yzt_fft(:,:,i)=fftshift(fft2(Ez_SI_record_position3_yzt(:,:,i),N_fft_x,N_fft_y));
    
    %Ex_sca_SI_record_yzt_fft(:,:,i)=fftshift(fft2(Ex_SI_record_yzt_sca(:,:,i),N_fft_x,N_fft_y));
    %Ey_sca_SI_record_yzt_fft(:,:,i)=fftshift(fft2(Ey_SI_record_yzt_sca(:,:,i),N_fft_x,N_fft_y));
    %Ez_sca_SI_record_yzt_fft(:,:,i)=fftshift(fft2(Ez_SI_record_yzt_sca(:,:,i),N_fft_x,N_fft_y));
    
    %Ex_inc_SI_record_yzt_fft(:,:,i)=fftshift(fft2(Ex_SI_record_yzt_vaccum(:,:,i),N_fft_x,N_fft_y));
    %Ey_inc_SI_record_yzt_fft(:,:,i)=fftshift(fft2(Ey_SI_record_yzt_vaccum(:,:,i),N_fft_x,N_fft_y));
    %Ez_inc_SI_record_yzt_fft(:,:,i)=fftshift(fft2(Ez_SI_record_yzt_vaccum(:,:,i),N_fft_x,N_fft_y));
end


for i=1:timesteps

    figure(1);
    subplot(2,3,1);
    mesh(double(Ez_SI_record_position1_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    title(['E,position1,t=',num2str(i)]);
    subplot(2,3,2);
    %mesh(double(Ez_SI_record_position2_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    pcolor(double(Ez_SI_record_position2_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    %caxis([0 1e10]);
    title(['E,position2,t=',num2str(i)]);
    subplot(2,3,3);
    mesh(double(Ez_SI_record_position3_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    title(['E,position3,t=',num2str(i)]);
    
    subplot(2,3,4);
    pcolor(F_y,F_x,(abs(double(Ez_SI_record_position1_yzt_fft(:,:,i)))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)');
    title(['fft_{E,position1},t=',num2str(i)]);
    temp=Ez_SI_record_position1_yzt_fft(:,:,i);
    caxis([0 color_factor*max(max(abs(double((temp(colormap_F))))))]);
    ylabel('k_y (k_0)');xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    subplot(2,3,5);
    pcolor(F_y,F_x,(abs(double(Ez_SI_record_position2_yzt_fft(:,:,i)))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)');
    title(['fft_{E,position2},t=',num2str(i)]);
    temp=Ez_SI_record_position2_yzt_fft(:,:,i);
    caxis([0 color_factor*max(max(abs(double((temp(colormap_F))))))]);
    ylabel('k_y (k_0)');xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    subplot(2,3,6);
    pcolor(F_y,F_x,(abs(double(Ez_SI_record_position3_yzt_fft(:,:,i)))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)');
    title(['fft_{E,position3},t=',num2str(i)]);
    temp=Ez_SI_record_position3_yzt_fft(:,:,i);
    caxis([0 color_factor*max(max(abs(double((temp(colormap_F))))))]);
    ylabel('k_y (k_0)');xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    k=waitforbuttonpress;    
    %figure(2);
    %subplot(2,3,1);
    %pcolor(double(theta_record_position1_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    %title(['theta,position1,t=',num2str(i)]); colorbar;
    %subplot(2,3,2);
    %pcolor(double(theta_record_position2_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    %title(['theta,position2,t=',num2str(i)]); colorbar;
    %subplot(2,3,3);
    %pcolor(double(theta_record_position3_yzt(:,:,i)));shading interp;colormap(jet(2^10));
    %title(['theta,position3,t=',num2str(i)]); colorbar;
    
end

clear;clc;close all;
working_dictionary='C:/Data/Theses_Postdocs/Stage_ludovic/Projet_RI_Chopineau/Code-FDTD'
exci='/save_temp/';
load_pulse_num=1; %%% which pulse to load
total_pulse_num=1; % 1 or 2 or more
plot_load_pulse_num=0; %%% plot the load_pulse_num or not
all_arrange=[1,2]; %%% how to arrange the subplot
rough_thick_multi_mode=importdata([working_dictionary,exci,'rough_thick_multi_mode.mat']);
%rough_thick_multi_mode=0;

substracting_source_flag=0; %%% substracting the fourier components of the source, 1-yes, 0-no.

y_corr=importdata([working_dictionary,exci,'y_corr.mat']);
z_corr=importdata([working_dictionary,exci,'z_corr.mat']);
N_points_PML=importdata([working_dictionary,exci,'N_points_PML.mat']);
N_points_z_all=importdata([working_dictionary,exci,'N_points_z_all.mat']);
N_points_y_all=importdata([working_dictionary,exci,'N_points_y_all.mat']);
N_points_x_film=importdata([working_dictionary,exci,'N_points_x_film.mat']);
N_points_y_main=importdata([working_dictionary,exci,'N_points_y_main.mat']);
N_points_z_main=importdata([working_dictionary,exci,'N_points_z_main.mat']);
delta_x=importdata([working_dictionary,exci,'delta_x.mat']); %%% in meter
rough_thick=importdata([working_dictionary,exci,'rough_thick.mat']);

points_no_roughness_edge=importdata([working_dictionary,exci,'points_no_roughness_edge.mat']);
sor_air_edge=importdata([working_dictionary,exci,'sor_air_edge.mat']);
source_edge=importdata([working_dictionary,exci,'source_edge.mat']);
roughness_function=importdata([working_dictionary,exci,'roughness_function','_pulse_',num2str(load_pulse_num),'.mat']); %%% roughness function

rough_height=zeros(N_points_y_main+1,N_points_z_main+1); %%% acting as the function of an AFM
for i=1:N_points_y_main+1
    for j=1:N_points_z_main+1
        for k=N_points_x_film+1:-1:1
            if roughness_function(k,i,j)==1
                rough_height(i,j)=(N_points_x_film+1-k)*delta_x;
                %rough_height(i,j)=k*delta_x;
                break;
            end
        end
    end
end
%crop_factor=0.2;
%rough_height=rough_height(round(crop_factor*(N_points_y_main+1)):end-round(crop_factor*(N_points_y_main+1)),round(crop_factor*(N_points_z_main+1)):end-round(crop_factor*(N_points_z_main+1)));
Y_dimension_main=(N_points_y_main+1)*delta_x;
D_focal=(Y_dimension_main-2*source_edge*delta_x)/2;
gau_ord=5; %%% super gaussian distribution of order gau_ord, for gau_ord=2, it is a gaussian distribution.
w0_guess=D_focal/(4^(1/gau_ord)); %%% beam waist at the focus
YZ_rough=meshgrid(y_corr(N_points_PML+1+1:N_points_y_all-(N_points_PML+1)),z_corr(N_points_PML+1+1:N_points_z_all-(N_points_PML+1)));
y_shift_guess=y_corr(round(length(y_corr)/2));
z_shift_guess=z_corr(round(length(z_corr)/2));
window_function=exp(-(2*(abs((YZ_rough-y_shift_guess)).^gau_ord+abs((YZ_rough'-z_shift_guess)).^gau_ord))./(w0_guess^gau_ord));
rough_height=rough_height.*window_function; %%% multiply by a window function to increase the quality of the fft map

if load_pulse_num>=2
    for i=N_points_x_film+1:-1:1
        if min(min(roughness_function(i,:,:)))==1
            where_to_look=N_points_x_film+1-i-rough_thick-rough_thick_multi_mode; %%% search for the layer that is immediately below the new roughness
            break;
        else
            where_to_look=N_points_x_film-rough_thick-rough_thick_multi_mode; %%% where to look at in the depth, +0 is immediately below the roughness.
        end
    end
else
    where_to_look=2;
end
see_depth=rough_thick+rough_thick_multi_mode+where_to_look;

k_range=5; %%% k0, the range of wave number to see

margins_=[-2 -3 -0.5 0.1];papersize_=[18 18/2^0.5];font_size_=20;marker_size_=6;line_width_=0.5;

if plot_load_pulse_num==1
    
    E_ab=importdata([working_dictionary,exci,'E_ab','_pulse_',num2str(load_pulse_num),'.mat']); %%% absorbed energy in the film
    E_ab=double(E_ab); %%% in J/m^3
    Ex_ab=importdata([working_dictionary,exci,'Ex_ab','_pulse_',num2str(load_pulse_num),'.mat']); %%% absorbed x energy in the film
    Ex_ab=double(Ex_ab); %%% in J/m^3
    Ey_ab=importdata([working_dictionary,exci,'Ey_ab','_pulse_',num2str(load_pulse_num),'.mat']); %%% absorbed y energy in the film
    Ey_ab=double(Ey_ab); %%% in J/m^3
    Ez_ab=importdata([working_dictionary,exci,'Ez_ab','_pulse_',num2str(load_pulse_num),'.mat']); %%% absorbed z energy in the film
    Ez_ab=double(Ez_ab); %%% in J/m^3
    E_xplusy_ab=Ex_ab+Ey_ab; %%% in J/m^3
    
    Ex_0_film=importdata([working_dictionary,exci,'Ex_0_film','_pulse_',num2str(load_pulse_num),'.mat']); %%% V/m field amplitude of the last half optical cycle
    Ex_0_film=double(Ex_0_film);
    Ey_0_film=importdata([working_dictionary,exci,'Ey_0_film','_pulse_',num2str(load_pulse_num),'.mat']); %%% V/m field amplitude of the last half optical cycle
    Ey_0_film=double(Ey_0_film);
    Ez_0_film=importdata([working_dictionary,exci,'Ez_0_film','_pulse_',num2str(load_pulse_num),'.mat']); %%% V/m field amplitude of the last half optical cycle
    Ez_0_film=double(Ez_0_film);
    E_0_film=importdata([working_dictionary,exci,'E_0_film','_pulse_',num2str(load_pulse_num),'.mat']); %%% V/m field amplitude of the last half optical cycle
    E_0_film=double(E_0_film);
    
    Ez_ab=1e-6*double(gather(Ez_ab)); %%% in J/cm^3
    Ez_ab_slice=reshape(Ez_ab(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    E_ab=1e-6*double(gather(E_ab)); %%% in J/cm^3
    E_ab_slice=reshape(E_ab(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    Ex_ab=1e-6*double(gather(Ex_ab)); %%% in J/cm^3
    Ex_ab_slice=reshape(Ex_ab(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    Ey_ab=1e-6*double(gather(Ey_ab)); %%% in J/cm^3
    Ey_ab_slice=reshape(Ey_ab(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    E_xplusy_ab=1e-6*double(E_xplusy_ab); %%% in J/cm^3
    E_xplusy_ab_slice=reshape(E_xplusy_ab(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    
    Ex_0_film_slice=reshape(Ex_0_film(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    Ey_0_film_slice=reshape(Ey_0_film(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    Ez_0_film_slice=reshape(Ez_0_film(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    E_0_film_slice=reshape(E_0_film(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    
    
    roughness_function_first_layer=reshape(roughness_function(end,:,:),N_points_y_main+1,N_points_z_main+1);
    roughness_rough=roughness_function_first_layer(points_no_roughness_edge+1:end-points_no_roughness_edge-1,points_no_roughness_edge+1:end-points_no_roughness_edge-1); %%% only the rough part %%%
    size_rough=size(roughness_rough);
    
    if substracting_source_flag==1
        %%% get the fourier components of the soruce
        Y_dimension_main=(N_points_y_main+1)*delta_x;
        D_focal=(Y_dimension_main-2*source_edge*delta_x)/2;
        gau_ord=6; %%% super gaussian distribution of order gau_ord, for gau_ord=2, it is a gaussian distribution.
        w0_guess=D_focal/(4^(1/gau_ord)); %%% beam waist at the focus
        YZ=meshgrid(y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)));
        amp_guess_E_ab=max(E_ab_slice(:));
        amp_guess_Ez_ab=max(Ez_ab_slice(:));
        y_shift_guess=y_corr(round(length(y_corr)/2));
        z_shift_guess=z_corr(round(length(z_corr)/2));
        error_fun_E_ab=@(x0) x0(1)*exp(-(2*(abs((YZ-x0(2))).^gau_ord+abs((YZ'-x0(3))).^gau_ord))./(x0(4)^gau_ord))-E_ab_slice;
        error_fun_Ez_ab=@(x0) x0(1)*exp(-(2*(abs((YZ-x0(2))).^gau_ord+abs((YZ'-x0(3))).^gau_ord))./(x0(4)^gau_ord))-Ez_ab_slice;
        %[coeff_E_ab,~,~,~]=lsqnonlin(error_fun_E_ab,[amp_guess_E_ab y_shift_guess z_shift_guess w0_guess]);
        [coeff_Ez_ab,resnorm,residual,exitflag]=lsqnonlin(error_fun_Ez_ab,[amp_guess_Ez_ab y_shift_guess z_shift_guess w0_guess]);
        E_ab_slice_fitted=coeff_E_ab(1)*exp(-(2*(abs((YZ-coeff_E_ab(2))).^gau_ord+abs((YZ'-coeff_E_ab(3))).^gau_ord))./(coeff_E_ab(4)^gau_ord));
        Ez_ab_slice_fitted=coeff_Ez_ab(1)*exp(-(2*(abs((YZ-coeff_Ez_ab(2))).^gau_ord+abs((YZ'-coeff_Ez_ab(3))).^gau_ord))./(coeff_Ez_ab(4)^gau_ord));
        %%% get the fourier components of the soruce
    else
        E_ab_slice_fitted=0*E_ab_slice;
        Ez_ab_slice_fitted=0*Ez_ab_slice;
    end
    
end
%%% premeter for fft and ifft %%%
N_fft_x=max(2^10,N_points_y_main+1-2*source_edge);
N_fft_y=max(2^10,N_points_z_main+1-2*source_edge);

lambda=800e-9;
k_0=2*pi/lambda;
F_spatial=2*pi*(1/(delta_x))/k_0; %%% sampling frequency (spatial) normliazed to k_0
F_x=linspace(-F_spatial/2,F_spatial/2,N_fft_x); %%% coordinate in k space %%%
F_y=linspace(-F_spatial/2,F_spatial/2,N_fft_y); %%% coordinate in k space %%%
F_xy_mesh=meshgrid(F_x,F_y);

%%% use to defind the colormap for the fft image
F_colormap=1; %%% k_0
colormap_F=find(abs(F_xy_mesh)>F_colormap);
color_factor=0.5;
%%% use to defind the colormap for the fft image

if plot_load_pulse_num==1
    
    rough_fft=fftshift(fft2(roughness_rough,N_fft_x,N_fft_y)/(size_rough(1)*size_rough(2))); %%% fouirer transform of the roughness founction
    rough_height_fft=fftshift(fft2(rough_height,N_fft_x,N_fft_y));
    
    Ex_ab_slice_fft=fftshift(fft2(Ex_ab_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    yita_x=Ex_ab_slice_fft./abs(rough_fft); %%% efficacy factor?
    
    Ey_ab_slice_fft=fftshift(fft2(Ey_ab_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    yita_y=Ey_ab_slice_fft./abs(rough_fft); %%% efficacy factor?
    
    %Ez_ab_slice_fft=fftshift(fft2(Ez_ab_slice,N_fft_x,N_fft_y)/((N_points_y_main+1)*(N_points_z_main+1))); %%% fouirer transform of the E_ab_slice
    Ez_ab_slice_fft=fftshift(fft2(Ez_ab_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    yita_z=Ez_ab_slice_fft./abs(rough_fft); %%% efficacy factor?
    Ez_ab_slice_fitted_fft=fftshift(fft2(Ez_ab_slice_fitted,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice_fitted
    Ez_ab_slice_fft=Ez_ab_slice_fft-Ez_ab_slice_fitted_fft; %%% substracting the fourier components of the source
    
    E_ab_slice_fft=fftshift(fft2(E_ab_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    yita=E_ab_slice_fft./abs(rough_fft); %%% efficacy factor?
    E_ab_slice_fitted_fft=fftshift(fft2(E_ab_slice_fitted,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice_fitted
    E_ab_slice_fft=E_ab_slice_fft-E_ab_slice_fitted_fft; %%% substracting the fourier components of the source
    
    E_xplusy_ab_slice_fft=fftshift(fft2(E_xplusy_ab_slice,N_fft_x,N_fft_y));
    
    E_0_film_slice_fft=fftshift(fft2(E_0_film_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    Ex_0_film_slice_fft=fftshift(fft2(Ex_0_film_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    Ey_0_film_slice_fft=fftshift(fft2(Ey_0_film_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    Ez_0_film_slice_fft=fftshift(fft2(Ez_0_film_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    
    h_rough_fun=figure('name','roughness function');
    pcolor(roughness_rough);shading interp;colormap(jet(2^10));axis equal;
    
    h_rough_fft=figure('name','rough_fft');
    pcolor(F_y,F_x,(abs(double(rough_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 0.003]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_E_ab=figure('name','E_ab');
    pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),...
        1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),...
        double(E_ab_slice));
    shading interp;colormap(jet(2^10));axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_E_ab_fft=figure('name','fft_E_ab');
    pcolor(F_y,F_x,(abs(double(E_ab_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(E_ab_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_E_xplusy_ab_fft=figure('name','fft_E_xplusy_ab');
    pcolor(F_y,F_x,(abs(double(E_xplusy_ab_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(E_xplusy_ab_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_Ez_ab=figure('name','Ez_ab');
    pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),...
        1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),...
        double(Ez_ab_slice));
    shading interp;colormap(jet(2^10));axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_Ez_ab_fft=figure('name','fft_Ez_ab');
    pcolor(F_y,F_x,(abs(double(Ez_ab_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(Ez_ab_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_Ey_ab=figure('name','Ey_ab');
    pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),...
        1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),...
        double(Ey_ab_slice));
    shading interp;colormap(jet(2^10));axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_Ey_ab_fft=figure('name','fft_Ey_ab');
    pcolor(F_y,F_x,(abs(double(Ey_ab_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(Ey_ab_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_Ex_ab=figure('name','Ex_ab');
    pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),...
        1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),...
        double(Ex_ab_slice));
    shading interp;colormap(jet(2^10));axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_Ex_ab_fft=figure('name','fft_Ex_ab');
    pcolor(F_y,F_x,(abs(double(Ex_ab_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(Ex_ab_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_E_ab_slice_fitted_y_dirction=figure('name','E_ab_slice_fitted, y direction');
    plot(1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),E_ab_slice_fitted(:,round(length(z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)))/2))); hold on;
    plot(1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),E_ab_slice(:,round(length(z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)))/2)),'LineStyle','none','Marker','.');
    xlabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_E_0_film=figure('name','E_0_flim');
    pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),...
        1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),...
        E_0_film_slice);
    shading interp;colormap(jet(2^10));axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_Ez_0_film=figure('name','Ez_0_flim');
    pcolor(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),...
        1e6*y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)),...
        Ez_0_film_slice);
    shading interp;colormap(jet(2^10));axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    h_E_0_film_fft=figure('name','fft_E_0_film_fft');
    pcolor(F_y,F_x,(abs(double(E_0_film_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(E_0_film_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    hx_E_0_film_fft=figure('name','fft_Ex_0_film_fft');
    pcolor(F_y,F_x,(abs(double(Ex_0_film_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(Ex_0_film_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    hy_E_0_film_fft=figure('name','fft_Ey_0_film_fft');
    pcolor(F_y,F_x,(abs(double(Ey_0_film_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(Ey_0_film_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    hz_E_0_film_fft=figure('name','fft_Ez_0_film_fft');
    pcolor(F_y,F_x,(abs(double(Ez_0_film_slice_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(Ez_0_film_slice_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
    h_E_ab_slice_fitted_z_dirction=figure('name','E_ab_slice_fitted, z direction');
    plot(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),E_ab_slice_fitted(round(length(y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)))/2),:)); hold on;
    plot(1e6*z_corr(N_points_PML+1+1+source_edge:N_points_z_all-(N_points_PML+1+source_edge)),E_ab_slice(round(length(y_corr(N_points_PML+1+1+source_edge:N_points_y_all-(N_points_PML+1+source_edge)))/2),:),'LineStyle','none','Marker','.');
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    figure('name','roughness_function_yz');pcolor(1e6*z_corr(N_points_PML+1+1:N_points_z_all-(N_points_PML+1)),...
        1e6*y_corr(N_points_PML+1+1:N_points_y_all-(N_points_PML+1)),...
        reshape(roughness_function(end-see_depth,:,:),N_points_y_main+1,N_points_z_main+1));shading interp;axis equal;
    xlabel('z (\mu m)','Interpreter','tex','FontSize',font_size_);ylabel('y (\mu m)','Interpreter','tex','FontSize',font_size_);
    set(gca,'fontsize',font_size_,'FontAngle','italic');
    
    figure('name','roughness_function_xz');pcolor(reshape(roughness_function(:,round((N_points_y_main+1)/2),:),N_points_x_film+1,N_points_z_main+1));shading interp
    
    figure('name','roughness_function_xy');pcolor(reshape(roughness_function(:,:,round((N_points_z_main+1)/2)),N_points_x_film+1,N_points_y_main+1));shading interp
    
    figure('name','rough_height');pcolor(rough_height);shading interp
    
    h_rough_height_fft=figure('name','fft_rough_height');
    pcolor(F_y,F_x,(abs(double(rough_height_fft))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(rough_height_fft(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);axis equal;xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);
    
end

rough_height_fft_all=cell(1,total_pulse_num);
figure('name','rough_height_all_pulses');
for ii=1:total_pulse_num
    roughness_function=importdata([working_dictionary,exci,'roughness_function','_pulse_',num2str(ii),'.mat']); %%% roughness function
    rough_height=zeros(N_points_y_main+1,N_points_z_main+1); %%% acting as the function of an AFM
    for i=1:N_points_y_main+1
        for j=1:N_points_z_main+1
            for k=N_points_x_film+1:-1:1
                if roughness_function(k,i,j)==1
                    rough_height(i,j)=(N_points_x_film+1-k)*delta_x;
                    %rough_height(i,j)=k*delta_x;
                    break;
                end
            end
        end
    end
    Y_dimension_main=(N_points_y_main+1)*delta_x;
    D_focal=(Y_dimension_main-2*source_edge*delta_x)/2;
    gau_ord=5; %%% super gaussian distribution of order gau_ord, for gau_ord=2, it is a gaussian distribution.
    w0_guess=D_focal/(4^(1/gau_ord)); %%% beam waist at the focus
    YZ_rough=meshgrid(y_corr(N_points_PML+1+1:N_points_y_all-(N_points_PML+1)),z_corr(N_points_PML+1+1:N_points_z_all-(N_points_PML+1)));
    y_shift_guess=y_corr(round(length(y_corr)/2));
    z_shift_guess=z_corr(round(length(z_corr)/2));
    window_function=exp(-(2*(abs((YZ_rough-y_shift_guess)).^gau_ord+abs((YZ_rough'-z_shift_guess)).^gau_ord))./(w0_guess^gau_ord));
    rough_height=rough_height.*window_function; %%% multiply by a window function to increase the quality of the fft map
    rough_height_fft=fftshift(fft2(rough_height,N_fft_x,N_fft_y));
    rough_height_fft_all{ii}=rough_height_fft;
    subaxis(all_arrange(1),all_arrange(2),ii,'spacing',0,'Margin',0);pcolor(rough_height);shading interp;axis off;
end

E_ab_slice_fft_all=cell(1,total_pulse_num);
figure('name','E_ab_all_pulses');
for ii=1:total_pulse_num
    E_ab=importdata([working_dictionary,exci,'E_ab','_pulse_',num2str(ii),'.mat']); %%% absorbed energy in the film
    E_ab=1e-6*double(E_ab); %%% in J/cm^3
    if ii>=2
        for i=N_points_x_film+1:-1:1
            if min(min(roughness_function(i,:,:)))==1
                where_to_look=N_points_x_film+1-i-rough_thick-rough_thick_multi_mode; %%% search for the layer that is immediately below the new roughness
                break;
            else
                where_to_look=N_points_x_film-rough_thick-rough_thick_multi_mode; %%% where to look at in the depth, +0 is immediately below the roughness.
            end
        end
    else
        where_to_look=2;
    end
    see_depth=rough_thick+rough_thick_multi_mode+where_to_look;
    E_ab_slice=reshape(E_ab(end-see_depth,source_edge+1:end-source_edge,source_edge+1:end-source_edge),...
        N_points_y_main+1-2*source_edge,N_points_z_main+1-2*source_edge);
    E_ab_slice_fft=fftshift(fft2(E_ab_slice,N_fft_x,N_fft_y)); %%% fouirer transform of the E_ab_slice
    E_ab_slice_fft_all{ii}=E_ab_slice_fft;
    subaxis(all_arrange(1),all_arrange(2),ii,'spacing',0,'Margin',0);pcolor(E_ab_slice);shading interp;axis off;
end

save([working_dictionary,exci,'rough_height_fft_all.mat'],'rough_height_fft_all');
save([working_dictionary,exci,'E_ab_slice_fft_all.mat'],'E_ab_slice_fft_all');
save([working_dictionary,exci,'F_x.mat'],'F_x');
save([working_dictionary,exci,'F_y.mat'],'F_y');
figure('name','rough_height_fft_all_pulses');
for ii=1:total_pulse_num
    subaxis(all_arrange(1),all_arrange(2),ii,'spacing',0,'Margin',0);
    pcolor(F_y,F_x,(abs(double(rough_height_fft_all{ii}))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(rough_height_fft_all{ii}(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);axis off;
end

figure('name','E_ab_slice_fft_all_pulses');
for ii=1:total_pulse_num
    subaxis(all_arrange(1),all_arrange(2),ii,'spacing',0,'Margin',0);
    pcolor(F_y,F_x,(abs(double(E_ab_slice_fft_all{ii}))));shading interp;colormap(jet(2^10));xlabel('k_x (k_0)','Interpreter','tex','FontSize',font_size_);
    caxis([0 color_factor*max(max(abs(double(E_ab_slice_fft_all{ii}(colormap_F)))))]);
    ylabel('k_y (k_0)','Interpreter','tex','FontSize',font_size_);xlim([-k_range k_range]);ylim([-k_range k_range]);
    set(gca,'fontsize',font_size_,'FontAngle','italic','XTick',-k_range:1:k_range,'YTick',-k_range:1:k_range);axis off;
end

%fillPage(h1, 'margins',margins_, 'papersize', papersize_);
%saveas(h1,[working_dictionary,exci,'E_absor_z_yz'],'fig');
%print(h1,'-dpdf','-r400',[working_dictionary,exci,'E_absor_z_yz']);

%fillPage(h2, 'margins',margins_, 'papersize', papersize_);
%saveas(h2,[working_dictionary,exci,'fft_roughness'],'fig');
%print(h2,'-dpdf','-r400',[working_dictionary,exci,'fft_roughness']);

%fillPage(h3, 'margins',margins_, 'papersize', papersize_);
%saveas(h3,[working_dictionary,exci,'fft_E_absor_z_yz'],'fig');
%print(h3,'-dpdf','-r400',[working_dictionary,exci,'fft_E_absor_z_yz']);
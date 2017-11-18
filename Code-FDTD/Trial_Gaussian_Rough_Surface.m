%% Surface spectral density
%%
% W(k) is the average distribution of each spatial frequency component of
% the randomly fluctuated surface profile.

clear all; close all; clc;

lambda = 800; % en nm
h = lambda; % RMS height
l = lambda; % correlation length

k = linspace(-0.005, 0.005); % wave number

% For a Gaussian correlation, the corresponding spectral density is :
W_g = (h^2 * l * exp(-(k * l).^2/4)) / sqrt(4*pi);

% For a exponential correlation, the corresponding spectral density is :
W_e = (h^2 * l) ./ (sqrt(4*pi) * (1 + (k * l).^2));

figure;
plot(k, W_g, 'r'); title('Surface spectral density - Gaussian/Exponential'); xlabel('Wavenumber');
ylabel('abs(W)');
hold on;
grid on;
plot(k, W_e, 'b'); legend('Gaussian', 'Exponential');
hold off;

% Gaussian profil evolution in function of RMS height
% couleurs = hsv(5);
% figure;
% for i = 1:5
%     hold on;
%     grid on;
%     plot(k,((i*h)^2*l*exp(-(k*l).^2/4))/sqrt(4*pi), 'color', couleurs(i,:));
%     title('Gaussian profil evolution in function of RMS height');
%     xlabel('Wavenumber');
%     ylabel('abs(W)');
%     legend(num2str(i));
%     leg(i)=cellstr(num2str(i));
% end
% legend(leg)
%
% % Exponential profil evolution in function of RMS height
% figure;
% for i = 1:5
%     hold on;
%     grid on;
%     plot(k,((i*h)^2*l)./(sqrt(4*pi)*(1 + (k*l).^2)), 'color', couleurs(i,:));
%     title('Exponential profil evolution in fonction of RMS height');
%     xlabel('Wavenumber');
%     ylabel('abs(W)');
%     legend(num2str(i));
%     leg(i)=cellstr(num2str(i));
% end
% legend(leg)

% Gaussian profil evolution in function of correlation length
couleurs = hsv(5);
figure;
for i = 1:5
    hold on;
    grid on;
    plot(k, (h^2 * (i * l) * exp(-(k * (i * l)).^2/4))/sqrt(4*pi), 'color', couleurs(i, :));
    title('Gaussian profil evolution in function of correlation length');
    xlabel('Wavenumber');
    ylabel('abs(W)');
    legend(num2str(i));
    leg(i) = cellstr(['lx = ', num2str(i), ' lambda']);
end
legend(leg)

% Exponential profil evolution in function of correlation length
% figure;
% for i = 1:5
%     hold on;
%     grid on;
%     plot(k,(h^2*(i*l))./(sqrt(4*pi)*(1 + (k*(i*l)).^2)), 'color', couleurs(i,:));
%     title('Exponential profil evolution in fonction of correlation length');
%     xlabel('Wavenumber');
%     ylabel('abs(W)');
%     legend(num2str(i));
%     leg(i)=cellstr(num2str(i));
% end
% legend(leg)
%% Generation of rough surface 1D numerically

couleurs = hsv(3);
figure;
% lambda en nm
for i = 1:3
    hold on;
    grid on;
    [f, x] = Roughness_Defined_1D(2000, 50, 10, 5);
    plot(x, f, 'color', couleurs(i, :));
    title('1D Gaussian random rough surface with Gaussian autocovariance function');
    xlabel('Wavenumber');
    ylabel('height');
    legend(num2str(i));
    leg(i) = cellstr(num2str(i));
end
legend(leg)

figure;
for i = 1:3
    hold on;
    grid on;
    [f, x] = Roughness_Defined_1D(2000, 50, 10*i, 5);
    plot(x, f, 'color', couleurs(i, :));
    title('1D Gaussian random rough surface in function of RMS height');
    xlabel('Wavenumber');
    ylabel('height');
    legend(num2str(i));
    leg(i) = cellstr(num2str(i));
end
legend(leg)

figure;
for i = 1:3
    hold on;
    grid on;
    [f, x] = Roughness_Defined_1D(2000, 50, 10, 5*i);
    plot(x, f, 'color', couleurs(i, :));
    title('1D Gaussian random rough surface in function of correlation length');
    xlabel('Wavenumber');
    ylabel('height');
    legend(num2str(i));
    leg(i) = cellstr(num2str(i));
end
legend(leg)

%% Generation of rough surface 2D numerically

% Cas isotrope lx = ly
figure;
%for i = 1:6
grid on;
%[Z,x,y] = Roughness_Defined_2D(120,40,5,2,1); % en nm
[f, x, y] = Discretised_roughness(311, 4, 0.2/100, 6, 60);

%length(f)

%[f,x,y] = Roughness_Defined_2D(120,40,5,2,2);
%
New2 = [];
New1 = [];
%New0 = [];
for i = 1:311
    for j = 1:311
        if f(i, j) == 2
            New2(i, j) = 1;
            New1(i, j) = 1;
            %New0(i,j) = 1;
        elseif f(i, j) == 1
            New2(i, j) = 0;
            New1(i, j) = 1;
            %New0(i,j) = 1;
        else
            New2(i, j) = 0;
            New1(i, j) = 0;
            %New0(i,j) = 1;
        end
    end
end
New2;
New1;


%   subplot(2,3,i);
%plot(f);
title('2D Gaussian random rough surface anisotrope');
xlabel('X direction');
ylabel('Y direction');
zlabel('height');
%end

M = [];
for j = 1:120
    M(j) = f(35, j);
end
plot(M)
title('Coupe transversale à y = 35 pour ly = 10*lx');
xlabel('X');
ylabel('Amplitude');
%zlabel('height');
%%
figure;
for i = 1:6
    grid on;
    [f, x, y] = Roughness_Defined_2D(30, 3, 10*i, 5, 5);
    subplot(2, 3, i);
    surfc(x, y, f);
    title(['2D Gaussian random rough surface isotrope for h = ', num2str(i*10)]);
    xlabel('Y direction');
    ylabel('X direction');
    zlabel('height');
end
%%

figure;
for i = 1:6
    hold on;
    [f, x, y] = Roughness_Defined_2D(30, 3, 10, 5*i, 5*i);
    subplot(2, 3, i);
    surfc(x, y, f);
    title('2D Gaussian random rough surface isotrope in function of correlation length');
    xlabel('X direction');
    ylabel('Y direction');
    zlabel('height');
end

%%
% Cas anisotrope lx ~= ly
figure;
for i = 1:6
    grid on;
    [f, x, y] = Roughness_Defined_2D(30, 3, 10, 6*i, 5);
    subplot(2, 3, i);
    surfc(x, y, f)
    title(['2D Gaussian random rough surface anisotrope in function of lx = ', num2str(i*6)]);
    xlabel('Y direction');
    ylabel('X direction');
    zlabel('height');
end
%%

figure;
for i = 1:6
    grid on;
    [f, x, y] = Roughness_Defined_2D(30, 3, 10, 5, 6*i);
    subplot(2, 3, i);
    surfc(x, y, f);
    title('2D Gaussian random rough surface anisotrope in function of ly');
    xlabel('X direction');
    ylabel('Y direction');
    zlabel('height');
end
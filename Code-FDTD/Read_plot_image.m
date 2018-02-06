function [ I ] = Read_plot_image( computation_box_size_y,computation_box_size_z )
% Read_plot_image Utility function to read an image into a convenient
% format
% Retrieve an image
% Crop unecessary informations via an interactive tool
% Resize the image to computable image in case if it is too big

clc; close all;
[filename, pathname] = uigetfile( ...
    {'*.*', 'All Files (*.*)'}, ...
    'Pick the image to analyze');% Retrieve the filename and its path
MEB_name = strcat(pathname, filename);
[II, ~] = imread(MEB_name);
[~, ~, Iz] = size(II);
if Iz > 3
    Ip = II(:, :, 4); % Retrieve the depth map (?)
    II = II(:, :, 1:3);
end

% Remove the informations at bottom if necessary
display(sprintf('Select the relevant part of the image'));
Icut = imcrop(II);
if Iz >= 3 % non B&W image corner-case
    J = rgb2gray(im2double(Icut));
else
    J = im2double(Icut);
end
I = J(end:-1:1, :); % normalize

% Reconstruction of real spatial scales
[PixelWidth, PixelHeight] = getpixelsize(MEB_name,I);
I = checksize(I,computation_box_size_y,computation_box_size_z,PixelWidth, PixelHeight);

[lig, col] = size(I);
xx = linspace(-(col / 2)*PixelWidth, (col / 2)*PixelWidth, col);
yy = linspace(-(lig / 2)*PixelHeight, (lig / 2)*PixelHeight, lig);
[X, Y] = meshgrid(xx*1.e6, yy*1.e6); % For graphics in �m

% Display the cropped pictures
Imlect = figure('Name', 'Image read', 'NumberTitle', 'off');
pcolor(X, Y, I);
colormap parula; shading interp; colorbar
get(gca);
xlabel('X [�m]', 'FontSize', 20, 'fontweight', 'bold');
ylabel('Y [�m]', 'FontSize', 20, 'Rotation', 90, 'fontweight', 'bold');
set(gca, 'TickDir', 'out', 'XTick', [-10, -5, 0, 5, 10], 'YTick', [-10, -5, 0, 5, 10]);
set(gca, 'linewidth', 2);
bar = colorbar;
set(bar, 'linewidth', 2);
caxis([0.3, 0.6])
% Prog to be developed

[filename, pathname] = uigetfile( ...
    {'*.*', 'All Files (*.*)'}, ...
    'Pick a file');% Retrieve the filename and its path
MEB_name = strcat(pathname, filename);
[II, map] = imread(MEB_name);
[Ix, Iy, Iz] = size(II);
if Iz > 3
    Ip = II(:, :, 4); % Retrieve the depth map (?)
    II = II(:, :, 1:3);
end

% Remove the informations at bottom if necessary
Icut = imcrop(II);
% Center the image % if Iy>Ix
[Ix, Iy, Iz] = size(Icut);
xymin = min(Ix, Iy);
dif = (Iy - Ix);
if dif <= 0 % Square image corner-case
    dif = 1;
end
Ic = Icut(1:xymin, ceil(dif/2):ceil((Iy - dif)/2)-1, 1:Iz);
if Iz >= 3 % non B&W image corner-case
    J = rgb2gray(im2double(Ic));
end
I = J(end:-1:1, :) / max(max(J)); % normalize

% Reconstruction des �chelles spatiales
[PixelWidth, PixelHeight] = getpixelsize(MEB_name);
if (PixelWidth == 0) || (PixelHeight == 0)
    prompt = {'Enter the pixel width (in \mu m) :', 'Enter the pixel height (in \mu m):'};
    dlg_title = 'Enter the pixel heights';
    num_lines = 1;
    defaultans = {'1', '1'};
    options.Interpreter = 'tex';
    answer = inputdlg(prompt, dlg_title, num_lines, defaultans, options);
    PixelWidth = str2double(answer{1});
    PixelHeight = str2double(answer{2});
end
[lig, col] = size(I);
xx = linspace(-(col / 2)*PixelWidth, (col / 2)*PixelWidth, col);
yy = linspace(-(lig / 2)*PixelHeight, (lig / 2)*PixelHeight, lig);
[X, Y] = meshgrid(xx*1.e6, yy*1.e6); % For graphics in �m


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
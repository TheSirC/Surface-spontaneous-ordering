% Prog to be developed

Filepath_read=XXX;
MEB_name=strcat(XXX,suffixfile);
II=imread(MEB_name);
[Ix Iy Iz]=size(II);

%Remove the informations at bottom if necessary
Icut=II(1:Ix-nbcellinfo,1:Iy,1:Iz);
%Center the image % if Iy>Ix
[Ix Iy Iz]=size(Icut);
xymin=min(Ix,Iy);
dif=(Iy-Ix);
Ic = Icut(1:xymin,floor(dif/2):floor(Iy-dif/2)-1,1:Iz);
J = rgb2gray(im2double(Ic));
I=J(end:-1:1,:)/max(max(J)); %normalize

% Reconstruction des échelles spatiales
[PixelWidth,PixelHeight]=getpixelsize(MEB_name);
[lig,col]=size(I);
xx=linspace(-(col/2)*PixelWidth,(col/2)*PixelWidth,col);
yy=linspace(-(lig/2)*PixelHeight,(lig/2)*PixelHeight,lig);
[X,Y]=meshgrid(xx*1.e6,yy*1.e6);     % For graphics in µm


Imlect = figure (1);
pcolor(X,Y,I);
colormap(map1) ;shading interp; colorbar
get(gca);
xlabel( 'X [µm]', 'FontSize', 20,'fontweight' ,'bold');
ylabel( 'Y [µm]', 'FontSize', 20, 'Rotation', 90 ,'fontweight' ,'bold');
set(gca, 'TickDir', 'out', 'XTick', [-10 -5 0 5 10], 'YTick',[-10 -5 0 5 10]);
set(gca,'linewidth',2);
bar=colorbar;
set(bar,'linewidth',2);
caxis([0.3 0.6])
function [PixelWidth,PixelHeight] = getlengthpixels( image )
%getlengthpixels Ask for the user to select a set of pixels that he knows
%the length. 
% CAUTION : We suppose at this point that the pixels are square !
%   image is the image to process

cdata = image;
[rows, columns] = size(cdata);
xdata = [1, columns];
ydata = [1, rows];

[x, y] = ginput(2);
if (x(2) < x(1))
    [x(1), x(2)] = deal(x(2), x(1));
end
pixel_length = x(2) - x(1);

prompt = {'Enter the length of the segment in \mu m:'};
dlg_title = 'Enter the pixel dimensions';
num_lines = 1;
defaultans = {'1'};
options.Interpreter = 'tex';
options.Resize = 'on';
options.WindowStyle = 'normal';
answer = inputdlg(prompt, dlg_title, num_lines, defaultans, options);
PixelWidth = str2double(answer{1})*1e-6*pixel_length;
PixelHeight = str2double(answer{1})*1e-6*pixel_length;

end


function [PixelWidth, PixelHeight] = getpixelsize(filename, image)
%getpixelsize  returns pixel width and height of MEB images.
%   Images must be in TIFF format
%   Open an image to retrieve its pixelsize if embedded in the file else
%   returns 0
% Usage: [PixelWidth,PixelHeight] = getpixelsize(filename)
info = imfinfo(filename);
if exist('info.BitsPerSample', 'var') == 0 % B&W image corner-case
    nbBits = 8;
    StripOffset = 8;
else
    nbBits = sum(info.BitsPerSample);
    StripOffset = info.StripOffsets(1);
end
offset = info.Width * info.Height * nbBits / 8;
fid = fopen(filename, 'r');
fseek(fid, offset+StripOffset, 'bof'); % on se d�place apr�s l'image
while ~feof(fid)
    txt = fgetl(fid);
    if size(sscanf(char(txt), 'PixelWidth=%g'))
        PixelWidth = sscanf(txt, 'PixelWidth=%g');
    end
    if size(sscanf(char(txt), 'PixelHeight=%g'))
        PixelHeight = sscanf(txt, 'PixelHeight=%g');
    end
end

fclose(fid);

% Handling corner-case where the variables where not found
if exist('PixelWidth', 'var') == 0
    PixelWidth = 0;
end

if exist('PixelHeight', 'var') == 0
    PixelHeight = 0;
end

% Ask the user for the information
if (PixelWidth == 0) || (PixelHeight == 0)
    f=helpdlg('Enter the pixel dimensions (leave zero if you want to choose a known interval)');
    waitfor(f);
    prompt = {'Enter the pixel width (in \mu m) :', 'Enter the pixel height (in \mu m):'};
	dlg_title = 'Pixels dimensions';
    num_lines = 1;
    defaultans = {'0', '0'};
    options.Interpreter = 'tex';
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    answer = inputdlg(prompt,dlg_title, num_lines, defaultans, options);
    PixelWidth = str2double(answer{1})*1e-6;
    PixelHeight = str2double(answer{2})*1e-6;
end

% Fallback if the pixels width and height are still unknown
if (PixelWidth == 0) || (PixelHeight == 0)
    f=helpdlg('Select the pixels of a known length ', 'Select a known length');
    waitfor(f);
    [PixelWidth,PixelHeight] = getlengthpixels(image);
end

end
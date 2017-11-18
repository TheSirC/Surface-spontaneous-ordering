function [PixelWidth, PixelHeight] = getpixelsize(filename)
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
fseek(fid, offset+StripOffset, 'bof'); % on se déplace après l'image
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
end
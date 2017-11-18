%  GETPIXELSIZE  returns pixel width and height of MEB images.
% Images must be in TIFF format
%
% Usage: [PixelWidth,PixelHeight] = getpixelsize(filename)

function [PixelWidth, PixelHeight] = getpixelsize(filename)

info = imfinfo(filename);
offset = info.Width * info.Height * sum(info.BitsPerSample) / 8;
fid = fopen(filename, 'r');
fseek(fid, offset+info.StripOffsets(1), 'bof'); % on se déplace après l'image
while ~feof(fid)
    txt = fgetl(fid);
    if size(sscanf(txt, 'PixelWidth=%g'))
        PixelWidth = sscanf(txt, 'PixelWidth=%g');
    end
    if size(sscanf(txt, 'PixelHeight=%g'))
        PixelHeight = sscanf(txt, 'PixelHeight=%g');
    end
end
fclose(fid);

end
%% Copyright (C) 2014-2015 Carnë Draug <carandraug@octave.org>
%%
%% This program is free software; you can redistribute it and/or modify it under
%% the terms of the GNU General Public License as published by the Free Software
%% Foundation; either version 3 of the License, or (at your option) any later
%% version.
%%
%% This program is distributed in the hope that it will be useful, but WITHOUT
%% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
%% details.
%%
%% You should have received a copy of the GNU General Public License along with
%% this program; if not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn  {Function File} {} imcrop ()
%% @deftypefnx {Function File} {} imcrop (@var{img})
%% @deftypefnx {Function File} {} imcrop (@var{ind}, @var{cmap})
%% @deftypefnx {Function File} {} imcrop (@var{h})
%% @deftypefnx {Function File} {} imcrop (@dots{}, @var{rect})
%% @deftypefnx {Function File} {[@var{cropped}] =} imcrop (@dots{})
%% @deftypefnx {Function File} {[@var{cropped}, @var{rect}] =} imcrop (@dots{})
%% @deftypefnx {Function File} {[@var{x}, @var{y}, @var{cropped}, @var{rect}] =} imcrop (@dots{})
%% Crop image.
%%
%% Displays the image @var{img} in a figure window and waits for the user to
%% select two points defining a bounding box.  For an indexed image, a
%% corresponding colormap can be specified in @var{cmap}.  For multi-dimensional
%% images (each 2D image is concatenated in the 4th dimension), only the
%% first image is displayed.
%%
%% If no image data is given, uses the current figure or figure from graphics
%% handle @var{h}.
%%
%% Non-interactive usage is supported if the last input argument is 4 element
%% vector @var{rect} defining the region of interest.  The first two elements
%% specify the initial @var{x_ini} and @var{y_ini} coordinates, and the last
%% two the @var{width} and @var{height}, i.e.,
%% @code{@var{rect} = [@var{x_ini} @var{y_ini} @var{width} @var{height}]}.
%% Note how this the opposite of the majority of Octave indexing rules where
%% rows come before columns.
%%
%% Returns the @var{cropped} image and a vector @var{rect} with the
%% coordinates and size for @var{cropped}.  If more than 3 output arguments
%% are requested, also returns the @var{x} and @var{y} data that define
%% the coordinate system.
%%
%% @emph{Note}: the values in @var{rect} are not necessarily integer values
%% and can't always be used directly as index values for other images.  To
%% crop the same region from a multiple images of the same size, either using
%% a multi-dimensional image:
%%
%% @example
%% @group
%% nd_img = cat (4, img1, img2, img3, img4);
%% cropped = imcrop (nd_img);
%% cropped_1 = cropped(:,:,:,1);
%% cropped_2 = cropped(:,:,:,2);
%% cropped_3 = cropped(:,:,:,3);
%% cropped_4 = cropped(:,:,:,4);
%% @end group
%% @end example
%%
%% or multiple calls to @code{imcrop}:
%%
%% @example
%% @group
%% [cropped_1, rect] = imcrop (img1);
%% cropped_2         = imcrop (img2, rect);
%% cropped_3         = imcrop (img3, rect);
%% cropped_4         = imcrop (img4, rect);
%% @end group
%% @end example
%%
%% @seealso{impixel, imshow}
%% @end deftypefn

%% TODO not yet implemented
%% @deftypefnx {Function File} {} imcrop (@var{xData}, @var{yData}, @dots{})

function varargout = imcrop(varargin)

%% Screw Matlab and their over complicated API's~ How can we properly
%% parse all the possible alternative calls? See
%% http://www.youtube.com/watch?v=1oZWacjmYm8 to understand how such
%% API's develop.

%% There is no check for this things, anything is valid. We (Octave)
%% are at least checking the number of elements otherwise the input
%% parsing would be horrible.
valid_rect = @(x) numel(x) == 4;
valid_system = @(x) numel(x) == 2;

rect = [];
interactive = true; % is interactive usage
alt_system = false; % was an alternative coordinate system requested?
from_fig = false; % do we have the image data or need to fetch from figure?

if (nargin > 5)
    print_usage();
end

rect = [];
if (numel(varargin) > 1 && valid_rect(varargin{end}))
    interactive = false;
    rect = varargin{end};
    varargin(end) = [];
end

xdata = [];
ydata = [];
if (numel(varargin) > 2 && valid_system(varargin{1}) && valid_system(varargin{2}))
    %% requested messy stuff
    %% we should probably reuse part of what impixel does
    alt_system = true;
    xdata = varargin{1};
    ydata = varargin{2};
    varargin([1, 2]) = [];
    error('imcrop: messing around with coordinate system is not implemented');
end

%% After we remove all that extra stuff, we are left with the image
fnargin = numel(varargin);
if (fnargin > 2)
    print_usage();
elseif (fnargin == 0)
    %% use current figure
    from_fig = true;
    h = gcf();
    %% We check isscalar() because ishandle() accepts arrays of handles, and we
    %% check '~= 0' because 0 is and handle for the 'root figure' which is
    %% invalid for imcrop (see bug %42714).
elseif (fnargin == 1 && isscalar(varargin{1}) ...
        && varargin{1} ~= 0 && ishandle(varargin{1}))
    %% use specified figure
    from_fig = true;
    h = varargin{1};
elseif (interactive)
    %% leave input check to imshow
    h = nd_imshow(varargin{:});
elseif (isimage(varargin{1}))
    %% we have the image data and it's not interactive, so there is
    %% nothing to do. We only check the minimum in the image.
else
    print_usage();
end

if (from_fig)
    hax = get(h, 'currentaxes');
    himage = findobj(hax, 'type', 'image');
    if (~isempty(himage)) ...
            himage = himage(1);
    else
        error('imcrop: expect the current axes to contain an image')
    end
    cdata = get(himage, 'cdata');
    xdata = get(himage, 'xdata');
    ydata = get(himage, 'ydata');
else
    cdata = varargin{1};
    if (~alt_system)
        [rows, columns] = size(cdata);
        xdata = [1, columns];
        ydata = [1, rows];
    end
end

%% Finally, crop the image
if (interactive)
    [x, y] = ginput(2);
    if (x(2) < x(1))
        [x(1), x(2)] = deal(x(2), x(1));
    end
    if (y(2) < y(1))
        [y(1), y(2)] = deal(y(2), y(1));
    end
    rect = [x(1), y(1), x(2) - x(1), y(2) - y(1)];
end
i_ini = max(round([rect(1), rect(2)]), [1, 1]);
tmp = size(cdata);
i_end = min(round([rect(1) + rect(3), rect(2) + rect(4)]), tmp(1:2));
img = cdata(i_ini(2):i_end(2), i_ini(1):i_end(1), :, :); % don't forget RGB and ND images

%% Even the API for the output is complicated
if (nargout == 0 && interactive)
    figure();
    %% In case we have a colormap or something like that, use
    %% it again when displaying the cropped image.
    nd_imshow(img, varargin{2:end});
elseif (nargout < 3)
    varargout{1} = img;
    varargout{2} = rect;
else
    varargout{1} = xdata;
    varargout{2} = ydata;
    varargout{3} = img;
    varargout{4} = rect;
end

end

%% shadows core function to support ND image.  If we have one, use
%% the first frame only
function h = nd_imshow(varargin)
size(varargin{1});
h = imshow(varargin{1}(:, :, :, 1), varargin{2:end});
end

%% test typical non-interactive usage, grayscale image
%~test
%~ a = randi (255, [100 100]);
%~ rect = [20 30 3 5];
%~ assert (nthargout ([1 2], @imcrop, a, rect), {a(30:35, 20:23) rect});
%~ assert (nthargout (2, @imcrop, a, rect), rect);
%~ assert (nthargout ([3 4], 4, @imcrop, a, rect), {a(30:35, 20:23) rect});

%% test typical non-interactive usage, RGB image
%~test
%~ rgb = randi (255, [100 100 3]);
%~ rect = [20 30 3 5];
%~ assert (nthargout ([1 2], @imcrop, rgb, rect), {rgb(30:35, 20:23,:) rect});
%~ assert (nthargout (2, @imcrop, rgb, rect), rect);
%~ assert (nthargout ([3 4], 4, @imcrop, rgb, rect), {rgb(30:35, 20:23,:) rect});

%% test typical non-interactive usage, indexed image
%~test
%~ a = randi (255, [100 100]);
%~ rect = [20 30 3 5];
%~ cmap = jet (255);
%~ assert (nthargout ([1 2], @imcrop, a, cmap, rect), {a(30:35, 20:23) rect});
%~ assert (nthargout (2, @imcrop, a, cmap, rect), rect);
%~ assert (nthargout ([3 4], 4, @imcrop, a, cmap, rect), {a(30:35, 20:23) rect});

%% test typical non-interactive usage, logical image
%~test
%~ a = rand (100) > 0.5;
%~ rect = [20 30 3 5];
%~ assert (nthargout ([1 2], @imcrop, a, rect), {a(30:35, 20:23) rect});
%~ assert (nthargout (2, @imcrop, a, rect), rect);
%~ assert (nthargout ([3 4], 4, @imcrop, a, rect), {a(30:35, 20:23) rect});

%% 0 is the root figure (always true figure handle), so make sure we use
%% scalar 0 as image data, not as figure handle.
%~assert (imcrop (0, [0.5 0.5 0.9 0.9]), 0);
%~assert (imcrop (zeros (5), [1 1 1 1]), zeros (2));

%% test out of bounds region of interest to crop (bug %49456)
%~test
%~ im = magic (5);
%~ assert (imcrop (im, [1 1 5 5]), im)
%~ assert (imcrop (im, [0 0 5 5]), im)
%~ assert (imcrop (im, [1 1 2 5]), im(:,1:3))
%~ assert (imcrop (im, [1 -3 2 5]), im(1:2,1:3))
%~ assert (imcrop (im, [5 -3 2 5]), im(1:2,5))

%~test
%~ %% Matlab returns [] (size 0x0) for this cases, while we return
%~ %% [] (size 2x0).  We are not compatible by design.  If it ever
%~ %% becomes an issue to anyone we can review this decision.
%~ assert (imcrop (magic (5), [6 -3 2 5]), zeros (2, 0))

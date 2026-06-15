%% SundewStemTracker.m
% reads in video showing the dynamics of sundew stem morphology changes
% during capture and converts them into a set of x,y point uniformly spaced
% along the stem vs time suitable for further analysis

%% initialization
close all; clear; clc;  % cleans up old windows, workspace, command widow

%  select which filter to use for the image
color_channel = 'Y';            % green; 'Y' = R+G = yellow

write_video = 1;
npts = 50;                      % number of points along stem backbone after resampling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the video file

[file path] = uigetfile('*.*','Select the experimental video file (.mp4 or .avi)');

videofile = [path file]
[video_pathstr,video_name,video_ext] = fileparts(videofile);

% create a video reader
% read the video frame

videoInReader = VideoReader(videofile);

% get total number of frames in the original video
nframes = videoInReader.NumberOfFrames;
framerate = videoInReader.FrameRate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the first and last frame to save as image files
prompt = {'First frame to save:','Last frame to save:','delta frame:', ...
    'Crop? 1=Y/0=N' 'framerate(fps)' 'pixel/mm' '# points to resample along stem'};
dlg_title = 'Video Analysis setups';
num_lines = 1;
def = {'1',num2str(nframes),'1','1',num2str(framerate),'1',num2str(npts)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
% Use curly bracket for subscript
first_frame = str2num(answer{1});
last_frame  = str2num(answer{2});
delta_frame = str2num(answer{3});
cropit      = str2num(answer{4});
framerate   = str2num(answer{5});
pixelpermm  = str2num(answer{6});
npts        = str2num(answer{7});

%% Construct a VideoWriter object, which creates a Motion-JPEG AVI file by default.
if write_video
    output_video_name = fullfile(video_pathstr,[video_name '_tracked.avi']);
    % create a movie object
    outputVideo = VideoWriter(output_video_name); 
    outputVideo.FrameRate = framerate;
    open(outputVideo);
end

%%
% loop through the video reading the images every delta_frame frames from
% first to last frame to analyze
% create figure and tiled layout for image and instructions
ij = 0;
for ii=first_frame:delta_frame:last_frame

    ij = ij + 1;                    % index for each tracked frame
    t(ij) = ii*framerate;           % time in s

    %% read in video frame by frame
    original_image = read(videoInReader,ii);
    % reduce haze
    original_image = imreducehaze(original_image);
    
    %% not sure how to use the different color channels
    %  The yellow stem is enhanced in the green 2nd channel of the color 
    %  RGB image (hairs are reddish, drops are white, background 
    %  brownish-grey)
    %  Subtracting the blue channel also removes background, 
    %  mount and drops nicely
    img = original_image(:,:,2) - original_image(:,:,3);

    %% crop if desired to region of interest

    croprgbimg = original_image;
    if ii==first_frame
        if cropit
            % crop the image and get its cropping box rect to use later
            [~,rect] = imcrop(original_image);
            % the crop rectangle parameters are: [xmin ymin width height]
            rect = round(rect);
            croprgbimg = imcrop(original_image,rect); close;
            nrows = rect(4); ncols = rect(3);   % use for correcting y origin later
            rgbimg = croprgbimg;   % save for nice plotting later
        else
            [nrows,ncols,ncolors] = size(img);  % use for correcting y origin later
        end
        % create meshgrid the same size as the image for mapping coordinagtes
        [X,Y] = meshgrid(1:ncols,1:nrows); 
    else
        if cropit
            croprgbimg = imcrop(original_image,rect);
            rgbimg = croprgbimg;   % for nice plotting later
        end
    end


    %% use the different color channels to get only the stem & enhance 
    % contrast
    %  The stem is yellow, Y = R + G, the hairs are reddish, drops white, 
    %  background brownish-grey), so using either G or Y works the best for
    %  discrimating the stem from other features in the  image

    switch color_channel
        case 'Y'
            imgRG = croprgbimg; imgRG(:,:,3) = 0;
            img = rgb2gray(imgRG);  % create a yellow = R + G image
        case 'G'
            img = croprgbimg(:,:,2);
    end

    % switch color_channel
    %     case 'R'
    %         img = croprgbimg(:,:,1);   % red channel only
    % 
    % end
    
   %% contrast enhancement 
    % to enhance contrast of the original image adapthisteq works much
    % better than histeq and slightly better than imadjust 
    img = adapthisteq(img,"NumTiles",[8 8]);

    %% select points along backbone
    figure('Position',[10 50 1800 500])
    tiledlayout(1, 2);
    % print instructions
    nexttile(1);
    axis off; % Hide axis lines for text area
    text(0.5, 0.5, {'Goal: choose enough points to define the stem backbone' ...
        'Click the left mouse button (normal button clicks) to add points', ...
        'A marker will appear where you click' ...
        'Remove the previously selected point by pressing the Backspace or Delete key' ...
        'Finish the selection using one of the following methods:'...
        'Double-click the left mouse button' ...
        'Shift-click or right-click' 'Press the Return or Enter key'}, ...
        'FontSize', 12);
    % show image and get points
    nexttile(2);
    imshow(img);
    [xi,yi] = getpts;
    % Click the left mouse button (normal button clicks) to add points. A marker will appear where you click.
    % Remove the previously selected point by pressing the Backspace or Delete key.
    % Finish the selection using one of the following methods:
    % 
    % Double-click the left mouse button.
    % Shift-click or right-click.
    % Press the Return or Enter key (finishes without adding a final point). 

    %% associate each point in xi,yi with a point on the meshgrid:
    xygrid = matchxytomeshgrid(xi,yi,X,Y);

    %% resample to get equally space points along the stem
    
    % first we need to remove any replicated points because
    % the code is fussy about successive replicate points:
    % to prevent errors, add tiny step between those with zero
    % difference before computing
    xyz_diff = diff(xygrid(:,1:2));  % get difference between points along backbone
    % get indices where both x and y are repeated
    idx_repl = find(xyz_diff(:,1) == 0 & xyz_diff(:,2) == 0);
    % add tiny random number to one of the identical
    % neighboring coordinates to distinguish repeated values
    xygrid(idx_repl,1:2) = xygrid(idx_repl,1:2) + ...
        1e-12*rand(numel(idx_repl),2);

    Pspline = interparc(npts,xygrid(:,1),xygrid(:,2),'spline');
    % compute evenly sampled points in units of mm:
    x_s = Pspline(:,1)/pixelpermm; 
    y_s = Pspline(:,2);
    close all;
    imshow(img); hold on; plot(x_s,y_s);plot(xi,yi,'o');
    axis equal; 
    % correct y for convention whereby the origin for y is at top left:
    y_s = (nrows - Pspline(:,2))/pixelpermm;
    
    %% save in the stem backbone array

    xy_stem(ij,:) = [t(ij) reshape([x_s y_s].', 1, [])];

    if write_video
        F = getframe(gcf);      % get the frame
        % write the frame to video
        writeVideo(outputVideo,F);%add the frames to the avi object created previously
    end
    close all;
end
display('done with reading in images, cropping & saving them');

outputfilename = fullfile(video_pathstr,[video_name '_stemxy.csv']);
header = {'time(s)'};
for i = 1:npts
    header = [header ['pt',num2str(i),'x(m)'] ['pt',num2str(i),'y(m)']];
end

outputarray = [header;num2cell(xy_stem)];
writecell(outputarray,outputfilename);

% close output video
if write_video, close(outputVideo), end

% end of main program

%% ************************************************************************
%% supporting functions
%% *************************************************************************
function closest_point = matchxytomeshgrid(xi,yi,X,Y)
    % for the set of path points xi yi,find the closest match in 
    % meshgrid XY = closest_point(1), closest_point(2)
    
    % based on Google Gemini code
    % Reshape the grid coordinates into a single list of points
    grid_points = [X(:), Y(:)];
    
    % Define the target points
    target_point = [xi yi];
    
    % Use knnsearch to find the closest grid point, ignoring distance
    [idx,~] = knnsearch(grid_points, target_point);
    
    % Get the coordinates of the closest point
    closest_point = grid_points(idx, :);
    
    % Display the result
    % fprintf('The closest grid point to (%f, %f) is (%f, %f) with a distance of %f\n', ...
    %     target_point(1), target_point(2), closest_point(1), closest_point(2), dist);
end
%% ************************************************************************
%% ************************************************************************
function [arclen,seglen] = arclength(px,py,pz,varargin)
% arclength: compute arc length of a space curve, or any curve represented as a sequence of points
% usage: [arclen,seglen] = arclength(px,py)         % a 2-d curve
% usage: [arclen,seglen] = arclength(px,py,pz)      % a 3-d space curve
% usage: [arclen,seglen] = arclength(px,py,method)  % specifies the method used
%
% Computes the arc length of a function or any
% general 2-d, 3-d or higher dimensional space
% curve using various methods.
%
% arguments: (input)
%  px, py, pz, ... - vectors of length n, defining points
%        along the curve. n must be at least 2. Replicate
%        points should not be present in the curve.
%
%  method - (OPTIONAL) string flag - denotes the method
%        used to compute the arc length of the curve.
%
%        method may be any of 'linear', 'spline', or 'pchip',
%        or any simple contraction thereof, such as 'lin',
%        'sp', or even 'p'.
%
%        method == 'linear' --> Uses a linear chordal
%               approximation to compute the arc length.
%               This method is the most efficient.
%
%        method == 'pchip' --> Fits a parametric pchip
%               approximation, then integrates the
%               segments numerically.
%
%        method == 'spline' --> Uses a parametric spline
%               approximation to fit the curves, then
%               integrates the segments numerically.
%               Generally for a smooth curve, this
%               method may be most accurate.
%
%        DEFAULT: 'linear'
%
%
% arguments: (output)
%  arclen - scalar total arclength of all curve segments
%
%  seglen - arclength of each independent curve segment
%           there will be n-1 segments for which the
%           arc length will be computed.
%
%
% Example:
% % Compute the length of the perimeter of a unit circle
% theta = linspace(0,2*pi,10);
% x = cos(theta);
% y = sin(theta);
%
% % The exact value is
% 2*pi
% % ans =
% %          6.28318530717959
%
% % linear chord lengths
% arclen = arclength(x,y,'l')
% % arclen =
% %           6.1564
%
% % Integrated pchip curve fit
% arclen = arclength(x,y,'p')
% % arclen =
% %          6.2782
%
% % Integrated spline fit
% arclen = arclength(x,y,'s')
% % arclen =
% %           6.2856
%
% Example:
% % A (linear) space curve in 5 dimensions
% x = 0:.25:1;
% y = x;
% z = x;
% u = x;
% v = x;
%
% % The length of this curve is simply sqrt(5)
% % since the "curve" is merely the diagonal of a
% % unit 5 dimensional hyper-cube.
% [arclen,seglen] = arclength(x,y,z,u,v,'l')
% % arclen =
% %           2.23606797749979
% % seglen =
% %         0.559016994374947
% %         0.559016994374947
% %         0.559016994374947
% %         0.559016994374947
%
%
% See also: interparc, spline, pchip, interp1
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/10/2010
% unpack the arguments and check for errors
if nargin < 2
    error('ARCLENGTH:insufficientarguments', ...
        'at least px and py must be supplied')
end

n = length(px);
% are px and py both vectors of the same length?
if ~isvector(px) || ~isvector(py) || (length(py) ~= n)
    error('ARCLENGTH:improperpxorpy', ...
        'px and py must be vectors of the same length')
elseif n < 2
    error('ARCLENGTH:improperpxorpy', ...
        'px and py must be vectors of length at least 2')
end
% compile the curve into one array
data = [px(:),py(:)];
% defaults for method and tol
method = 'linear';
% which other arguments are included in varargin?
if numel(varargin) > 0
    % at least one other argument was supplied
    for i = 1:numel(varargin)
        arg = varargin{i};
        if ischar(arg)
            % it must be the method
            validmethods = {'linear' 'pchip' 'spline'};
            ind = strmatch(lower(arg),validmethods);
            if isempty(ind) || (length(ind) > 1)
                error('ARCLENGTH:invalidmethod', ...
                    'Invalid method indicated. Only ''linear'',''pchip'',''spline'' allowed.')
            end
            method = validmethods{ind};

        else
            % it must be pz, defining a space curve in higher dimensions
            if numel(arg) ~= n
                error('ARCLENGTH:inconsistentpz', ...
                    'pz was supplied, but is inconsistent in size with px and py')
            end

            % expand the data array to be a 3-d space curve
            data = [data,arg(:)]; %#ok
        end
    end

end
% what dimension do we live in?
nd = size(data,2);
% compute the chordal linear arclengths
seglen = sqrt(sum(diff(data,[],1).^2,2));
arclen = sum(seglen);
% we can quit if the method was 'linear'.
if strcmpi(method,'linear')
    % we are now done. just exit
    return
end
% 'spline' or 'pchip' must have been indicated,
% so we will be doing an integration. Save the
% linear chord lengths for later use.
chordlen = seglen;
% compute the splines
spl = cell(1,nd);
spld = spl;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for i = 1:nd
    switch method
        case 'pchip'
            spl{i} = pchip([0;cumsum(chordlen)],data(:,i));
        case 'spline'
            spl{i} = spline([0;cumsum(chordlen)],data(:,i));
            nc = numel(spl{i}.coefs);
            if nc < 4
                % just pretend it has cubic segments
                spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
                spl{i}.order = 4;
            end
    end

    % and now differentiate them
    xp = spl{i};
    xp.coefs = xp.coefs*diffarray;
    xp.order = 3;
    spld{i} = xp;
end
% numerical integration along the curve
polyarray = zeros(nd,3);
for i = 1:spl{1}.pieces
    % extract polynomials for the derivatives
    for j = 1:nd
        polyarray(j,:) = spld{j}.coefs(i,:);
    end

    % integrate the arclength for the i'th segment
    % using quadgk for the integral. I could have
    % done this part with an ode solver too.
    seglen(i) = quadgk(@(t) segkernel(t),0,chordlen(i));
end
% and sum the segments
arclen = sum(seglen);
% ==========================
%   end main function
% ==========================
%   begin nested functions
% ==========================
    function val = segkernel(t)
        % sqrt((dx/dt)^2 + (dy/dt)^2)

        val = zeros(size(t));
        for k = 1:nd
            val = val + polyval(polyarray(k,:),t).^2;
        end
        val = sqrt(val);

    end % function segkernel
end % function arclength

%% ************************************************************************

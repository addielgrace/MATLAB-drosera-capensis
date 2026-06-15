%% initialization
close all; clear; clc;  % cleans up old windows, workspace, command widow
[file, path] = uigetfile('*.csv','Select the data file');
data_cellarray = readmatrix(fullfile(path, file));

%CHANGE THIS FOR EACH LEAF
dt_real = 95100 / (31); %in seconds
dt_realh = dt_real / 3600;
t_realh = (0:31)' * dt_realh;
mm_per_pixel = 34/ 353.9509;


npts_fit = 12;                           % number of points in the window fit to a circle
half_window = floor(npts_fit/2);

if mod(npts_fit,2) ~= 1                 % if it's even, add 1 -- npts_fit has to be odd
    npts_fit = npts_fit + 1;
end

% remove header
data_cellarray(1,:) = [];

%find x and y coordinates in dataset
xdata = data_cellarray(:, 2:2:end);
xdata(:,2:end) = xdata(:,2:end) - xdata(:,1) +  xdata(1,1); %from suzanne, align all data to the same reference point to plot
xdata(:,1) = xdata(1,1);
xdata = xdata * mm_per_pixel;

ydata = data_cellarray(:,3:2:end);
ydata(:,2:end) = ydata(:,2:end) - ydata(:,1) +  ydata(1,1);
ydata(:,1) = ydata(1,1);
ydata = ydata * mm_per_pixel;

nframes = size(xdata,1);
npts = size(xdata,2);

for i = 1:height(data_cellarray)
    dx = diff(xdata(i,:));
    dy = diff(ydata(i,:));
    ds = sqrt(dx.^2 + dy.^2);
    s_all{i} = [0, cumsum(ds)];
    s = s_all{1}; %refernce for acrlength is the base of the leaf (first coordinate)
    s_norm_all{i} = s/s(end);  % normalized arc length (0 → 1)
end


%normalized arclength
s_norm = s_norm_all{1};
s_norm;

%"real life" arch length
s_ref = s_all{1}; % arclength from first frame

% find center coordinates by using Suzanne's curvature code and
% circlefitbypratt
centerxy = nan(nframes, npts, 2);
for i = 1:nframes
    XY = [xdata(i,:)' ydata(i,:)'];  % create an [npts,2] XY array for each frame
    for j = (half_window+1):(npts-half_window)
        Par = CircleFitByPratt(XY(j-half_window:j+half_window,:));

        centerx(i,j) = Par(1);  % x-center
        centery(i,j) = Par(2);  % y-center

    end
end



%make sure centerx, centery, and the x and y coordinates all have the same
%length by erasing data of first 4 data points (related to half_window, how
%Suzanne's code computed center points)
xdata_angle = xdata(:,half_window+1:end); % x coordinates 
ydata_angle = ydata(:, half_window+1:end); % y coordinates

% %define variables
xdata_angle1 = xdata_angle(1:end-1,:);
xdata_angle2 = xdata_angle(2:end,:);

ydata_angle1 = ydata_angle(1:end-1,:);
ydata_angle2 = ydata_angle(2:end,:);

centerx1 = centerx(1:end-1,:);
centerx2 = centerx(2:end,:);
centery1 = centery(1:end-1,:);
centery2 = centery(2:end,:);

%find angles of corresponding coordinates between different time frames
x1 = xdata_angle1 - centerx1; 
y1 = ydata_angle1 - centery1;
x2 = xdata_angle2 - centerx2;
y2 = ydata_angle2 - centery2;
angle1 = atan2(y1, x1); %find angles corresponding to uneven time frames (with respect to x axis)
angle1 = unwrap(angle1, [], 2);
angle2 = atan2(y2, x2); %find angles corresponding to even time frames (with respect to x axis)
angle2 = unwrap(angle2,[], 2);
angle3 = angle2 - angle1; % dtheta, subtracting the 2 angles

% % find angular velocity
omega = (angle3 ./(dt_realh)) ;

% Find midpoint between frame1 and frame2, since angular velocity is also
% technically calculated for the midpoints (dtheta/dt)
x_mid = 0.5 * (xdata_angle1 + xdata_angle2);
y_mid = 0.5 * (ydata_angle1 + ydata_angle2);


% Only plot certain timeframes
% x_p   = x_mid(1:2:end, :);
% y_p   = y_mid(1:2:end, :);
% omega_p   = omega(1:end, :);

% %plot all timeframes 
x_p = x_mid;
y_p  = y_mid;
omega_p = omega;
t_mid = t_realh;

nPoints = size(omega_p, 2);

% Match s_norm to trimmed data (because of half_window)
s_used = s_norm(1:nPoints);

x_p1 = x_p(1:3:end, :);
y_p1 = y_p(1:3:end, :);
omega_p1 = omega_p(1:3:end, :);

% plot curve plot
% all curves show the velocity at the midpoint between 2 time frames, so it
% represents the velocity between the earlier frame and the later frame
figure (1); hold on; axis equal
set(gcf, 'Units', 'pixels', 'Position', [1, 1, 657, 573]);
xlabel("x", "FontSize", 15, 'FontName', 'Arial');
ylabel("y","FontSize", 15, 'FontName', 'Arial');
title("Angular Velocity Mapped Onto Leaf Curves", 'FontSize', 19, 'FontName', 'Arial');
cb = colorbar;
c.Label.String = "Angular Velocity (rad/hour)";
c.Label.FontSize = 14;
c.Label.FontName = 'Arial';


ax = gca;
ax.FontSize = 14;      % increases tick numbers size
ax.FontName = 'Arial'; % optional, for consistency


for j = 1:size(omega_p1,1)

    x = x_p1(j,:).';
    y = y_p1(j,:).';
    wj = omega_p1(j,:).';

    surface([x x], [y y], zeros(length(x),2), [wj wj], ...
        'EdgeColor','interp', ...
        'FaceColor','none', ...
        'LineWidth', 2);
end

%plot heatmap
figure(2)
set(gcf, 'Units', 'pixels', 'Position', [1, 1, 657, 573]);
imagesc(s_norm,t_mid, omega_p)
% imagesc(s_norm, t_mid, omega_p, 'AlphaData', ~isnan(omega_p))
xlabel('Normalized Arc Length', "FontSize", 15, 'FontName', 'Arial' );
ylabel('Time (hours)', "FontSize", 15, 'FontName', 'Arial');
title('Angular Velocity over Time and Arc Length', "FontSize", 19, "FontName", "Arial");
c = colorbar;
c.Label.String = "Angular Velocity (rad/hour)";
c.Label.FontSize = 15;
c.Label.FontName = 'Arial';

ax = gca;
ax.FontSize = 14;      % increases tick numbers size
ax.FontName = 'Arial'; % optional, for consistency
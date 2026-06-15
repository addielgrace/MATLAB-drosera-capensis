close all; clear; clc;
[file, path] = uigetfile('*.csv','Select the data file');
data = readmatrix(fullfile(path, file));

mm_per_pixel = 16/633.1840;
dt_real = 12300 / (5-1); %in seconds
dt_realh = dt_real / 3600;
t_realh = (0:4)' * dt_realh;

x_coords = data(1:end, 2:2:end);
x_coords = x_coords * mm_per_pixel;

%%FIND FRONT VELOCITY BY USING HIGHEST Y COORDINATE
% y_coords = data(1:end, 3:2:end);
% y_coords = y_coords * mm_per_pixel;
% 
% %find arclength
% for i = 1:size(x_coords,1)
% 
%     dx1 = diff(x_coords(i,:));
%     dy1 = diff(y_coords(i,:));
%     ds = sqrt(dx1.^2 + dy1.^2);
%     s = [0, cumsum(ds)]; 
%     s_all{i} = s;
%     s_norm_all{i} = s/s(end);  % normalized arc length (0 → 1)
% end
% 
% %normalized arclength
% s_norm = s_norm_all{1};
% 
% %"real life" arch length
% s_ref = s_all{1} ;  % arclength from first frame
% 
% %find index and highest value for each row in y_coords
% [val,idx] = max(y_coords, [],2);
% FRONT VELOCITY
% just select one point along the length of the leaf (highest one)
% use that to compute veloicty
%find velocity
% vx = diff(x_coords,1,1)/dt_realh;
% vx1 = vx(:,1);
% 
% vy = diff(y_coords,1,1)/dt_realh;
% 
% v_mag = sqrt(vx.^2 + vy.^2);

% FRONT VELOCITY ONE POINT
vx = diff(x_coords,1,1)/dt_realh;
vx1 = vx(:,1)

%%SAVE IF FRONT VELOCITY IS IN 2 PARTS
% save('25leaf1part1.mat', 'vx1');
% 
% vmean1 = mean(vx1(1:6))
% sd1 = std(vx1(1:6))
% vmean2 = mean(vx1(7:end))
% sd2 = std(vx1(7:end))
vmean3 = mean(vx1(1:end))
sd3 = std(vx1(1:end))


figure(1);
plot(vx1(:,1),t_realh(2:end));
xlabel('horizontal component of front velocity (vx, mm/hour)' )
ylabel('time (hours)')

% figure(2);
% plot(vx1(:,2),t_realh(2:end));
% xlabel('horizontal component of front velocity (vx, mm/hour)' )
% ylabel('time (hours)')


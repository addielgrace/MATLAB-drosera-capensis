close all; clear; clc;
[file, path] = uigetfile('*.csv','Select the data file');
data = readmatrix(fullfile(path, file));

%extract time column from data set
time = data(:, 1);

%CHANGE THIS FOR EACH LEAF
dt_real = 95100 / (31); %in seconds
dt_realh = dt_real / 3600;
t_realh = (0:31)' * dt_realh;
mm_per_pixel = 34/ 353.9509;

x_ref = data(1,2);   % base point reference
y_ref = data(1,3);


x_coords = data(:, 2:2:end);
x_coords(:,2:end) = x_coords(:,2:end) - x_coords(:,1) +  x_coords(1,1);
x_coords(:,1) = x_coords(1,1);
x_coords = x_coords * mm_per_pixel;

y_coords = data(:, 3:2:end);
y_coords(:,2:end) = y_coords(:,2:end) - y_coords(:,1) +  y_coords(1,1);
y_coords(:,1) = y_coords(1,1);
y_coords = y_coords * mm_per_pixel;


%find arclength
for i = 1:height(data)

    dx = diff(x_coords(i,:));
    dy = diff(y_coords(i,:));
    ds = sqrt(dx.^2 + dy.^2);
    s = [0, cumsum(ds)]; 
    s_all{i} = s;
    s_norm_all{i} = s/s(end);  % normalized arc length (0 → 1)
end

%normalized arclength
s_norm = s_norm_all{1}

%"real life" arch length
s_ref = s_all{1}   % arclength from first frame

%find velocity
dx = diff(x_coords,1,1)/dt_realh;
dy = diff(y_coords,1,1)/dt_realh;
size(dx)

v_mag = sqrt(dx.^2 + dy.^2);

v_space_mean = mean(v_mag, 2)
% 
% % %delete outlier (when leaf falls)
% size(v_space_mean)
% v_space_mean(7) = [];

v_space_mean;

%Total translational velocity
vmean = mean(v_space_mean)

%initial translational velocity (change 1:x, depends on leaf)
vmean1 = mean(v_space_mean(1:19))

%subsequent translational velocity (change x:2, depends on leaf)
vmean2 = mean(v_space_mean(20:end))


figure(1)
imagesc(s_norm, t_realh(2:end), v_mag)
set(gca,'YDir','normal')
xlabel('Arclength')
ylabel('Time (hours)')
title('Velocity Plot')
c = colorbar;
c.Label.String = 'velocity';

figure(2); hold on; axis equal
xlabel("x"); ylabel("y");
title("Translational Velocity Mapped Onto Leaf Curves - Leaf 2");

cb = colorbar;
ylabel(cb, "translational velocity");


for j = 1:size(v_mag,1)

    x = x_coords(j+1,:).'; % VELOCITY SHOWN CORRESPONDS TO THE LATER CURVE
    y = y_coords(j+1,:).';
    vj = v_mag(j,:).';

    surface([x x], [y y], zeros(length(x),2), [vj vj], ...
        'EdgeColor','interp', ...
        'FaceColor','none', ...
        'LineWidth', 2);
end



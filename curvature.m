close all; clear; clc;
[file, path] = uigetfile('*.csv','Select the data file');
data = readmatrix(fullfile(path, file));

%CHANGE THIS FOR EACH LEAF
dt_real = 95100 / (31); %in seconds
dt_realh = dt_real / 3600;
t_realh = (0:31)' * dt_realh;
mm_per_pixel = 34/ 353.9509;



x_coords = data(:, 2:2:end);
x_coords(:,2:end) = x_coords(:,2:end) - x_coords(:,1) +  x_coords(1,1);
x_coords(:,1) = x_coords(1,1);
x_coords = x_coords * mm_per_pixel;

y_coords = data(:, 3:2:end);
y_coords(:,2:end) = y_coords(:,2:end) - y_coords(:,1) +  y_coords(1,1);
y_coords(:,1) = y_coords(1,1);
y_coords = y_coords * mm_per_pixel;

nrows = size(x_coords,1); % after deleting row
nPts = size(x_coords,2);
k_mat = nan(nrows, nPts);

%find arclength
for i = 1:size(x_coords,1)

    dx = diff(x_coords(i,:));
    dy = diff(y_coords(i,:));
    ds = sqrt(dx.^2 + dy.^2);
    s = [0, cumsum(ds)]; 
    s_all{i} = s;
    s_norm_all{i} = s/s(end);  % normalized arc length (0 → 1)


    % find first and second derivatives
    dx1 = gradient(x_coords(i,:), s);
    dx2 = gradient(dx1, s);
    dy1 = gradient(y_coords(i,:), s);
    dy2 = gradient(dy1, s);

    % calculate and store curvature (k)
    k = abs(dx1.*dy2 - dx2 .* dy1)./(dx1.^2 + dy1.^2).^(3/2);

    k_mat(i,:) = k;
    %integrate
    L = s_all{i}(end);
    K_total(i) = trapz(s, k);   % integrated curvature
    K_norm(i) = K_total(i)/L

end

%normalized arclength
s_norm = s_norm_all{1};

%"real life" arch length
s_ref = s_all{1} % arclength from first frame

% changecurvature = diff(k_mat,1,1)/dt_realh;
% 
% mean_dkdt_all = mean(abs(changecurvature), "all") %mean magnitude of curvature change
% 
% mean_dkdt_initial = mean(abs(changecurvature(1:1:4,:)), "all")
% 
% mean_dkdt_subsequent = mean(abs(changecurvature(5:1:end,:)), "all")

figure (1);
set(gcf, 'Units', 'pixels', 'Position', [1, 1, 657, 573]);
imagesc(s_norm, t_realh, k_mat);
load('myCustomColormap.mat')
colormap(CustomColormap)
cb = colorbar;
cb.Label.String = 'Curvature (1/mm)';
cb.Label.FontSize = 23;
cb.Label.FontName = 'Arial';
clim([0 1.5])
set(gca, 'FontSize', 23, 'FontName', 'Arial')
set(gca, 'YDir', 'normal')
xlabel('Normalized Arc Length', 'FontSize', 23, 'FontName', 'Arial')
ylabel('Time (hours)', 'FontSize', 23, 'FontName', 'Arial')


% eliminate some times from the data set to produce clearer figure
% x_coordsp = x_coords(1:3:end, :);
% y_coordsp = y_coords(1:3:end, :);
% k_mat_plot = k_mat(1:3:end, :);

% don't eliminate soe times from the data set to produce clearer figure
x_coordsp = x_coords;
y_coordsp = y_coords;
k_mat_plot = k_mat;


%transpose matrices to produce correct figure
x_coordst = x_coordsp' ;
y_coordst = y_coordsp';

figure (2); hold on; axis equal
set(gcf, 'Units', 'pixels', 'Position', [1, 1, 657, 573]);
colormap(CustomColormap)
cb = colorbar;
cb.Label.String = 'Curvature (mm^-1)';
cb.Label.FontSize = 23;
cb.Label.FontName = 'Arial';
cb = colorbar;
clim([0 1.5]);
ylabel(cb, "Curvature (1/mm)");
xlabel('x', 'FontSize', 23, 'FontName', 'Arial'); 
ylabel('y', 'FontSize', 23, 'FontName', 'Arial');
set(gca, 'FontSize', 23, 'FontName', 'Arial')

% xlim([90 290]);
% ylim([70 260]);

for j = 1:size(x_coordsp,1)   % loop over the plotted time steps (every 3rd frame)

    x = x_coordsp(j,:).';
    y = y_coordsp(j,:).';
    k = k_mat_plot(j,:).';

    % colored curve (color = k)
    surface([x x], [y y], zeros(length(x),2), [k k], ...
        'EdgeColor','interp', 'FaceColor','none', 'LineWidth',2);

end


% figure(3)
% imagesc(s_norm, t_realh, changecurvature);
% set(gca,'YDir','normal');
% xlabel( 'arclength')
% ylabel('time (hours)')
% title("Curvature heatmap");
% 
% size(t_realh)
% size(K_norm)


figure(3)
set(gcf, 'Units', 'pixels', 'Position', [1, 1, 657, 573]);
plot(t_realh, K_norm.', 'LineWidth',2)
set(gca, 'FontSize', 22, 'FontName', 'Arial')
xlabel('Time (hours)', 'FontSize', 23, 'FontName', 'Arial')
ylabel('Integrated Curvature (1/mm)', 'FontSize', 23, 'FontName', 'Arial')

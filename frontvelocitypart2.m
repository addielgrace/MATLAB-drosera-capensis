close all; clear; clc;

%COMBINE FREE VELOCITY OF MULITPLE PARTS IF NEEDED

% Load your saved velocity files
load('25leaf1part1.mat'); % this loads vx1 (SAVED 1)
vx_dataset1 = vx1 

load('25leaf1part2.mat'); % replace with your second MAT file
vx_dataset2 = vx1  % make sure the variable name matches what is inside the MAT file (SAVED 2)

% load('75leaf4part3.mat'); % replace with your second MAT file
% vx_dataset3 = vx1;  % make sure the variable name matches what is inside
% the MAT file (SAVED 3, IF NEEDED)

% Combine the datasets vertically
vx_combined = [vx_dataset1; vx_dataset2]

%COMBINE NORMALIZATION CONSTANTS
dt_real =  33600 / (11); %in seconds
dt_realh = dt_real / 3600;
t_realh = (0:11)' * dt_realh;


vmean1 = mean(vx_combined(1:8))
sd1 = std(vx_combined(1:8))
vmean2 = mean(vx_combined(9:end))
sd2 = std(vx_combined(9:end))
vmean3 = mean(vx_combined(1:end))
sd3 = std(vx_combined(1:end))

% Plot
figure;
plot(vx_combined(:,1), t_realh(2:end)); % first column as front velocity
xlabel('Time (hours)');
ylabel('Horizontal velocity vx (mm/hour)');
title('Combined Front Velocity from Two Datasets');
grid on;

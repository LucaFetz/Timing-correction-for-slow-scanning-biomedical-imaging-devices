%% clear workspace and add library
clc
clear all
close all
path = genpath('GlobalBioIm');
addpath(path);

%% define parameters
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 1; % 0 or 1
param.noise_snr = 13; %signal to noise ratio
param.lambda =0.00464159;% 2*3e-2; %thune regulation term. 0 for no regulation.
param.opti_type = "FBS"; %choose optimization algorithm
param.regul_type = "L1"; %choose regulation type
param.plot_flag = 0; %vilsualize steps of algo
param.GTconfig = 4; %different ground truths
param.GT.speed = 0.7;

param.samples_coordinates_x = "classic"; %coordinates for x sampling
param.samples_coordinates_y = "classic"; %coordinates for y sampling


%% assess reconstruction effect with different dynamics.

%% ------ GT4 ------
%fix noise at a given level
param.noise = 1;
param.noise_snr = 16;
repeats = 5;

%choose best regul and opti type
param.opti_type = "FBS"; %choose optimization algorithm
param.regul_type = "L1"; %choose regulation type
result_GradDsct_L2 = Etape_2D(param);

%choose speed between 0 and pi and lambda in logspace
speed_list = 0:0.5:3;
lambda_list = logspace(-5,1,7);

snr_measurements_FBS_L1 = zeros(length(speed_list), 1);
snr_reconstruction_FBS_L1 = zeros(length(speed_list), 1);

%assess reconstruction efficiency with different speeds and lambda
for m = 1:length(speed_list)
    param.GT.speed = speed_list(m);
    for n = 1:length(lambda_list)
        param.lambda = lambda_list(n);
        for o = 1:repeats
            result_FBS_L1(m,n) = Etape_2D(param);
            snr_measurements_FBS_L1(m,n) = snr_measurements_FBS_L1(m,n) + result_FBS_L1(m,n).snr_measurements/repeats;
            snr_reconstruction_FBS_L1(m,n) = snr_reconstruction_FBS_L1(m,n) + result_FBS_L1(m,n).snr_reconstruction/repeats;
            fprintf('test %d of %d. \n', (m-1)*length(lambda_list)+n, length(speed_list)*length(lambda_list))
        end
    end
end

figure
imagesc(snr_reconstruction_FBS_L1);
xlabel('lambda');ylabel('speed');
xticklabels(lambda_list);yticklabels(speed_list);
title('SNR of reconstructed frames')

figure
imagesc(snr_measurements_FBS_L1);
xlabel('lambda');ylabel('speed');
xticklabels(lambda_list);yticklabels(speed_list);
title('SNR of measurements')
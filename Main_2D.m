%% clear workspace and add library
clc
clear all
close all
path = genpath('GlobalBioIm');
addpath(path);

%% define parameters
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 0; % 0 or 1
param.noise_snr = 10; %signal to noise ratio
param.lambda =0;% 2*3e-2; %thune regulation term. 0 for no regulation.
param.opti_type = "GradDsct"; %choose optimization algorithm
param.regul_type = "L2"; %choose regulation type
param.plot_flag = 0; %vilsualize steps of algo
param.GTconfig = 2; %different ground truths

param.samples_coordinates_x = "backnforth"; %coordinates for x sampling
param.samples_coordinates_y = "backnforth"; %coordinates for y sampling

%% run ETAPE_2D
result = Etape_2D(param);

%% assess noise and lambda effects
repeats = 1; %number of repeats to mean noise effects and have a more robust estimation
param.noise = 1;
lambda_list = logspace(-4,0,20);
snr_list = 1:1:10;
snr_measurements = zeros(length(lambda_list), length(snr_list));
snr_reconstruction = zeros(length(lambda_list), length(snr_list));

for m = 1:length(lambda_list)
    param.lambda = lambda_list(m);
    for n = 1:length(snr_list)
        param.noise_snr = snr_list(n);
        for o = 1:repeats
            result(m,n) = Etape_2D(param);
            snr_measurements(m,n) = snr_measurements(m,n) + result(m,n).snr_measurements/repeats;
            snr_reconstruction(m,n) = snr_reconstruction(m,n) + result(m,n).snr_reconstruction/repeats;
            fprintf('test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
        end
    end
end

figure
imagesc(snr_reconstruction);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of reconstructed frames')

figure
imagesc(snr_measurements);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of measurements')

%% visualize best results
[best_m_idx, best_n_idx] = find(snr_reconstruction == max(max(snr_reconstruction)));
param.noise_snr = snr_list(best_n_idx); 
param.lambda = lambda_list(best_m_idx); 
param.plot_flag = 1; 
Etape_2D(param);
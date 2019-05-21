%% clear workspace and add library
clc
clear all
close all
path = genpath('GlobalBioIm');
addpath(path);

%% define parameters
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 0; % 0 or 1
param.noise_snr = 16; %signal to noise ratio
param.lambda =0.00464159;% 2*3e-2; %thune regulation term. 0 for no regulation.
param.opti_type = "FBS"; %choose optimization algorithm
param.regul_type = "TV"; %choose regulation type
param.plot_flag = 0; %vilsualize steps of algo
param.GTconfig = 4; %different ground truths
param.GT.speed = 1.5;

param.samples_coordinates_x = "classic"; %coordinates for x sampling
param.samples_coordinates_y = "classic"; %coordinates for y sampling

%% run ETAPE_2D
result_GradDsct_L2 = Etape_2D(param);

%% assess noise and lambda effects
% param.opti_type = "GradDsct"; %choose optimization algorithm
% param.regul_type = "L2"; %choose regulation type
% 
% snr_measurements = zeros(length(lambda_list), length(snr_list));
% snr_reconstruction = zeros(length(lambda_list), length(snr_list));
% 
% for m = 1:length(lambda_list)
%     param.lambda = lambda_list(m);
%     for n = 1:length(snr_list)
%         param.noise_snr = snr_list(n);
%         for o = 1:repeats
%             result(m,n) = Etape_2D(param);
%             snr_measurements(m,n) = snr_measurements(m,n) + result(m,n).snr_measurements/repeats;
%             snr_reconstruction(m,n) = snr_reconstruction(m,n) + result(m,n).snr_reconstruction/repeats;
%             fprintf('test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
%         end
%     end
% end
% 
% figure
% imagesc(snr_reconstruction);
% xlabel('noise');ylabel('regul');
% xticklabels(snr_list);yticklabels(lambda_list);
% title('SNR of reconstructed frames')
% 
% figure
% imagesc(snr_measurements);
% xlabel('noise');ylabel('regul');
% xticklabels(snr_list);yticklabels(lambda_list);
% title('SNR of measurements')

%% assess noise and lambda effects with different optimization algos and regul types
repeats = 5; %number of repeats to mean noise effects and have a more robust estimation
param.noise = 1;
lambda_list = logspace(-5,1,7);
snr_list = 16:1:16;

%---- GradDsct and L2 -----
param.opti_type = "GradDsct"; %choose optimization algorithm
param.regul_type = "L2"; %choose regulation type

snr_measurements_GradDsct_L2 = zeros(length(lambda_list), length(snr_list));
snr_reconstruction_GradDsct_L2 = zeros(length(lambda_list), length(snr_list));

for m = 1:length(lambda_list)
    param.lambda = lambda_list(m);
    for n = 1:length(snr_list)
        param.noise_snr = snr_list(n);
        for o = 1:repeats
            result_GradDsct_L2(m,n) = Etape_2D(param);
            snr_measurements_GradDsct_L2(m,n) = snr_measurements_GradDsct_L2(m,n) + result_GradDsct_L2(m,n).snr_measurements/repeats;
            snr_reconstruction_GradDsct_L2(m,n) = snr_reconstruction_GradDsct_L2(m,n) + result_GradDsct_L2(m,n).snr_reconstruction/repeats;
            fprintf('L2 test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
        end
    end
end

figure
imagesc(snr_reconstruction_GradDsct_L2);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of reconstructed frames L2')

figure
imagesc(snr_measurements_GradDsct_L2);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of measurements')

%---- FBS and L1 -----
param.opti_type = "FBS"; %choose optimization algorithm
param.regul_type = "L1"; %choose regulation type

snr_measurements_FBS_L1 = zeros(length(lambda_list), length(snr_list));
snr_reconstruction_FBS_L1 = zeros(length(lambda_list), length(snr_list));

for m = 1:length(lambda_list)
    param.lambda = lambda_list(m);
    for n = 1:length(snr_list)
        param.noise_snr = snr_list(n);
        for o = 1:repeats
            result_FBS_L1(m,n) = Etape_2D(param);
            snr_measurements_FBS_L1(m,n) = snr_measurements_FBS_L1(m,n) + result_FBS_L1(m,n).snr_measurements/repeats;
            snr_reconstruction_FBS_L1(m,n) = snr_reconstruction_FBS_L1(m,n) + result_FBS_L1(m,n).snr_reconstruction/repeats;
            fprintf('L1 test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
        end
    end
end

figure
imagesc(snr_reconstruction_FBS_L1);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of reconstructed frames L1')

%---- FBS and TV -----
param.opti_type = "FBS"; %choose optimization algorithm
param.regul_type = "TV"; %choose regulation type

snr_measurements_FBS_TV = zeros(length(lambda_list), length(snr_list));
snr_reconstruction_FBS_TV = zeros(length(lambda_list), length(snr_list));

for m = 1:length(lambda_list)
    param.lambda = lambda_list(m);
    for n = 1:length(snr_list)
        param.noise_snr = snr_list(n);
        for o = 1:repeats
            result_FBS_TV(m,n) = Etape_2D(param);
            snr_measurements_FBS_TV(m,n) = snr_measurements_FBS_TV(m,n) + result_FBS_TV(m,n).snr_measurements/repeats;
            snr_reconstruction_FBS_TV(m,n) = snr_reconstruction_FBS_TV(m,n) + result_FBS_TV(m,n).snr_reconstruction/repeats;
            fprintf('TV test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
        end
    end
end

figure
imagesc(snr_reconstruction_FBS_TV);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of reconstructed frames TV')



%---- FBS and L12 -----
param.opti_type = "FBS"; %choose optimization algorithm
param.regul_type = "L12"; %choose regulation type

snr_measurements_FBS_L12 = zeros(length(lambda_list), length(snr_list));
snr_reconstruction_FBS_L12 = zeros(length(lambda_list), length(snr_list));

for m = 1:length(lambda_list)
    param.lambda = lambda_list(m);
    for n = 1:length(snr_list)
        param.noise_snr = snr_list(n);
        for o = 1:repeats
            result_FBS_L12(m,n) = Etape_2D(param);
            snr_measurements_FBS_L12(m,n) = snr_measurements_FBS_L12(m,n) + result_FBS_L12(m,n).snr_measurements/repeats;
            snr_reconstruction_FBS_L12(m,n) = snr_reconstruction_FBS_L12(m,n) + result_FBS_L12(m,n).snr_reconstruction/repeats;
            fprintf('L12 test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
        end
    end
end

figure
imagesc(snr_reconstruction_FBS_L12);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of reconstructed frames L12')

%---- Tikhonov and L2 -----
param.opti_type = "GradDsct"; %choose optimization algorithm
param.regul_type = "Tikhonov"; %choose regulation type

snr_measurements_GradDsct_Tikhonov = zeros(length(lambda_list), length(snr_list));
snr_reconstruction_GradDsct_Tikhonov = zeros(length(lambda_list), length(snr_list));

for m = 1:length(lambda_list)
    param.lambda = lambda_list(m);
    for n = 1:length(snr_list)
        param.noise_snr = snr_list(n);
        for o = 1:repeats
            result_GradDsct_Tikhonov(m,n) = Etape_2D(param);
            snr_measurements_GradDsct_Tikhonov(m,n) = snr_measurements_GradDsct_Tikhonov(m,n) + result_GradDsct_Tikhonov(m,n).snr_measurements/repeats;
            snr_reconstruction_GradDsct_Tikhonov(m,n) = snr_reconstruction_GradDsct_Tikhonov(m,n) + result_GradDsct_Tikhonov(m,n).snr_reconstruction/repeats;
            fprintf('Tikhonov test %d of %d. \n', (m-1)*length(snr_list)+n, length(lambda_list)*length(snr_list))
        end
    end
end

figure
imagesc(snr_reconstruction_GradDsct_Tikhonov);
xlabel('noise');ylabel('regul');
xticklabels(snr_list);yticklabels(lambda_list);
title('SNR of reconstructed frames Tikhonov')

%% visualize best results
% [best_m_idx, best_n_idx] = find(snr_reconstruction_GradDsct_L2 == max(max(snr_reconstruction_GradDsct_L2)));
% param.noise_snr = snr_list(best_n_idx); 
% param.lambda = lambda_list(best_m_idx); 
% param.plot_flag = 1; 
% Etape_2D(param);
%% clear workspace and add library
clc
clear all
close all
path = genpath('GlobalBioIm');
addpath(path);

%% define parameters
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 1; % 0 or 1
param.noise_snr = 10; %signal to noise ratio
param.lambda = 2*3e-2; %thune regulation term. 0 for no regulation.

%% Create 2D ground truth
%define GT parameters
param.GT.sigma_x = 2;
param.GT.centering_x = 3;
param.GT.sigma_y = 2;
param.GT.centering_y = 4;
param.GT.sigma_t = 2;
param.GT.centering_t = 5;
%create GT
[measurements, f0, GT, Nx, Ny, Nt] = create_2d_ground_truth(param.GT.sigma_x,param.GT.sigma_y,param.GT.sigma_t,param.GT.centering_x,param.GT.centering_y, param.GT.centering_t, param.noise, param.noise_snr);

%% plot GT with animation in time
figure;
for i = 1:Nt
    imagesc(GT(:,:,i));colorbar;caxis([0 1]);
    title('GT')
    pause(0.2);
end

%% plot measurements with animation in time
%since the reshape is column wise we reshape as [Ny,Nx,Nt] and then permute
measurements_matrix = permute(reshape(measurements,[Ny,Nx,Nt]),[2 1 3]);
figure;
for i = 1:Nt
    imagesc(measurements_matrix(:,:,i));colorbar;caxis([0 1]);
    title('measurements/naive approach')
    pause(0.2);
end

%% Create forward model
[H, h] = create_2d_forward_model(Nx,Ny,Nt,param.spline_order);

%% plot H
figure, imagesc(H);

%% Find C with inverse problem
C = optimize_c_2D(H, param.lambda, measurements');
C = permute(reshape(C',[Ny,Nx,Nt]),[2 1 3]);

%% plot C
figure
for i = 1:Nt
    imagesc(C(:,:,i));colorbar;caxis([0 1]);
    title('Coefs');
    pause(0.2);
end
%% Reconstruct good frames
[f,reconstructed_frames] = interpolate_2D(h, C);
%% visualize reconstructed frames
figure;
for i = 1:Nt
    imagesc(reconstructed_frames(:,:,i));colorbar;caxis([0 1]);
    title('reconstructed frames');
    pause(0.2);
end

%% assess reconstruction quality and visualize results
result.snr_measurements = snr(GT,GT-measurements_matrix);
result.snr_reconstruction = snr(GT,GT-reconstructed_frames);

figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:Nt

    str = sprintf('Time = %d', i)
    suptitle(str)
    
    subplot(231),imagesc(GT(:,:,i));colorbar;caxis([0 1]);
    title('GT')
    
    subplot(232),imagesc(measurements_matrix(:,:,i));colorbar;caxis([0 1]);
    title('measurements/naive approach')
    
    subplot(233),imagesc(reconstructed_frames(:,:,i));colorbar;caxis([0 1]);
    title('reconstructed frames');
    
    subplot(235),imagesc(GT(:,:,i) - measurements_matrix(:,:,i));colorbar;caxis([0 max(max(max(GT(:,:,:) - measurements_matrix(:,:,:))))]);
    title('GT - measurements');
    
    subplot(236),imagesc(GT(:,:,i) - reconstructed_frames(:,:,i));colorbar;caxis([0 max(max(max(GT(:,:,:) - reconstructed_frames(:,:,:))))]);
    title('GT - reconstructed frames');
    
    pause(0.2);
end

for i = 1:Nt

    str = sprintf('Time = %d', i)
    suptitle(str)
    
    subplot(231),surf(GT(:,:,i));colorbar;caxis([0 1]);zlim([0 1]);
    title('GT')
    
    subplot(232),surf(measurements_matrix(:,:,i));colorbar;caxis([0 1]);zlim([0 1]);
    title('measurements/naive approach')
    
    subplot(233),surf(reconstructed_frames(:,:,i));colorbar;caxis([0 1]);zlim([0 1]);
    title('reconstructed frames');
    
    subplot(235),surf(GT(:,:,i) - measurements_matrix(:,:,i));colorbar;caxis([min(min(min(GT(:,:,:) - measurements_matrix(:,:,:)))) max(max(max(GT(:,:,:) - measurements_matrix(:,:,:))))]);zlim([min(min(min(GT(:,:,:) - measurements_matrix(:,:,:)))) max(max(max(GT(:,:,:) - measurements_matrix(:,:,:))))]);
    title('GT - measurements');
    
    subplot(236),surf(GT(:,:,i) - reconstructed_frames(:,:,i));colorbar;caxis([min(min(min(GT(:,:,:) - reconstructed_frames(:,:,:)))) max(max(max(GT(:,:,:) - reconstructed_frames(:,:,:))))]);zlim([min(min(min(GT(:,:,:) - reconstructed_frames(:,:,:)))) max(max(max(GT(:,:,:) - reconstructed_frames(:,:,:))))]);
    title('GT - reconstructed frames');
    
    pause(0.2);
end
%% clear workspace and add library
clc
clear all
close all
path = genpath('GlobalBioIm-release');
addpath(path);

%% define parameters
spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
noise = 0; % 0 or 1
snr = 10; %signal to noise ratio
lambda = 0*3e-2; %thune regulation term. 0 for no regulation.

%% Create 2D ground truth
%define GT parameters
sigma_x = 2;
centering_x = 3;
sigma_y = 2;
centering_y = 4;
sigma_t = 2;
centering_t = 5;
%create GT
[measurements, f0, df0] = create_2d_ground_truth(sigma_x,sigma_y,sigma_t,centering_x,centering_y, centering_t);

%% Call correction function
%Frame_correction_2D(measurements,spline_order,noise,noise_snr,lambda)

%% visualize results
%% clear previous parameters and results
clear param result
clc
close all
path = genpath('GlobalBioIm');
addpath(path);

%% define parameters
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 1; % 0 or 1
param.snr = 10; %signal to noise ratio
param.lambda = 0*3e-2; %thune regulation term. 0 for no regulation.

%% run time correction once
Etape_1D_NoPlots;

%% execute time correction for different noise and regulation and store result
param.spline_order = 3;
param.noise = 1;
%faire plusieurs fois l'expérience
for m = 1:100
    param.lambda = m*3e-3; %utiliser logspace pour chercher lambda optimal
    for n = 1:10
        param.snr = 5 + n;
        Etape_1D_NoPlots;
        close all
        %result will have different Reg in each row and different noise in
        %each column
        result(m,n).df0 = df0;
        result(m,n).df = df;
        result(m,n).measurement = measurement;
        result(m,n).Fit = df0(:,1:Nx:Nx*Nt) - df(:,1:Nx:Nx*Nt);
        Error(m,n) = rms(result(m,n).Fit(:));
    end
end
%% look for best results for each noise value
[xTmp, xBest] = min(Error,[],1);
%% plot results
figure,imagesc(Error);

for i = 1:size(xBest,2)
    %comparison between sampled f, GT (at integer times) and measures in 3D
    figure, subplot(231),surf(result(xBest(i),i).df0(:,1:Nx:Nx*Nt));%zlim([0 1]);
    title('Sampled GT')
    xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
    subplot(232),surf(result(xBest(i),i).df(:,1:Nx:Nx*Nt));%zlim([0 1]);
    title('Reconstruction')
    xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
    subplot(233),surf(reshape(result(xBest(i),i).measurement,Nx,Nt));%zlim([0 1]);
    title('Measurements/Naive approach')
    xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
    subplot(235),surf(result(xBest(i),i).Fit);
    title('GT - Reconstruction')
    xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
    subplot(236),surf(result(xBest(i),i).df0(:,1:Nx:Nx*Nt) - reshape(result(xBest(i),i).measurement,Nx,Nt));
    title('GT - Measurements/Naive approach')
    xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
end
%%
tmpdf0 = result(xBest(10),10).df0(:,1:Nx:Nx*Nt);
tmpdf = result(xBest(10),10).df(:,1:Nx:Nx*Nt);
snr_reconstruction = 20*log10(norm(tmpdf0(:))/norm(tmpdf(:)-tmpdf0(:)));
%% 
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 1; % 0 or 1
param.snr = 15; %signal to noise ratio
param.lambda = 0*3e-2; %thune regulation term. 0 for no regulation.
Etape_1D_NoPlots;

%comparison between GT (at integer times), sampled f with and without
%Reg
figure, subplot(231),surf(result(xBest(10),10).df0(:,1:Nx:Nx*Nt));%zlim([0 1]);
title('Sampled GT')
xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
subplot(232),surf(result(xBest(10),10).df(:,1:Nx:Nx*Nt));%zlim([0 1]);
title('Reconstruction with Reg 0.1110')
xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
subplot(233),surf(df(:,1:Nx:Nx*Nt));%zlim([0 1]);
title('Reconstruction without Reg')
xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
subplot(235),surf(result(xBest(10),10).Fit);
title('GT - Reconstruction with Reg')
xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');
subplot(236),surf(result(xBest(10),10).df0(:,1:Nx:Nx*Nt) - df(:,1:Nx:Nx*Nt));
title('GT - Reconstruction without Reg')
xlabel('Time [frame]'), ylabel('Space [px]'), zlabel('Signal intensity');

tmpdf0 = df0(:,1:Nx:Nx*Nt);
tmpdf = df(:,1:Nx:Nx*Nt);
snr_reconstruction_noReg = 20*log10(norm(tmpdf0(:))/norm(tmpdf(:)-tmpdf0(:)));


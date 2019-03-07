%% clear workspace
clc
clear all
close all

%% Ground truth simulation
%continuous grounf truth f0
% define size and time ranges

% define variable and functions
f0 = @(x,t) gaussmf(x, [2 5]).*gaussmf(t, [2 5]); %continuous ground truth f0(x,t)

%Discrete ground truth df0
%define number of frames and size of image
Frames = 10;
Nt = Frames;
Nx = 11; %size in nb of pixels
sampling_time = 1/(Nx);
%sample continuous f0
dx = (0:1:Nx-1)'; 
%xImage = gaussmf(dx, [2 (Nx-1)/2])'; %1D image in x plane
%xImageStacked = repmat(xImage,1,Frames); %stack of Frames xImages
dt = (0:sampling_time:(Nx*Frames-1)*sampling_time)';
%t_fluctuation = gaussmf(dt, [2 5]); % fluctuation of xImage over time
%df0 = f0(repmat(dx',1,1),repmat(dt,1,Nx)); %sampled f0. each row is the image in a certain timeframe
df0 = f0(repmat(dx',Nx*Frames,Frames),repmat(dt,1,Nx*Frames)); %sampled f0. In reality only the first Nx columns are needed. This is used to simplify measurement with diag(f0)

%% measurement and usual approximation
measurement = diag(df0); %each sample is taken from its timeframe. In usual approximation it is considered as image at time t=0
figure
plot(repmat(dx,Frames,1), [measurement df0(1,:)' df0(floor(end/2),:)']);
legend(' superposition of Measurements / usual approximation','Ground truth at t=0','Ground truth at time t=1/2')
mean_measurement = mean(reshape(measurement,Nx,Frames),2);
figure
plot(dx, [mean_measurement df0(1,1:Nx)' df0(floor(end/2),1:Nx)']);
legend(' Mean of Measurements / usual approximation','Ground truth at t=0','Ground truth at time t=1/2')

figure 
for i=1:Frames
    plot(dx, [measurement((i-1)*Nx+1:i*Nx) df0(i,(i-1)*Nx+1:i*Nx)']);
    legend(' Measurement / usual approximation','Ground truth')
    title(i)
    pause(0.001)
end

%% forward model
%continuous forward model
k = repmat(0:1:Nx-1,Nt,1);
l = repmat((0:1:Nt-1)',1,Nx);
h = @(x,t) (B_spline(t-l).*B_spline(x-k));

for i = 0:1:Nx-1
    for j = 0:1:Nt-1 
       H(i*(Nx-1)+(j+1),:) = reshape(h(i,j+i*sampling_time),1,Nx*Nt); %does sampling_time have a meaning here?
    end
end
y = diag(measurement);
%% coefficient optimization for B spline interpolation
%c = inv(H).*y;
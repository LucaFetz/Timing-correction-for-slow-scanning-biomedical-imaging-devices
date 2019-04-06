function [measurements, f0, df0] = create_2d_ground_truth(sigma_x,sigma_y,sigma_t,centering_x,centering_y, centering_t)
%function that creates 2D+time ground truth(GT) as a 2D gaussian in space and a
%gaussian in time
%input: sigma and centering to define the gaussians
%output:measurements[(Nx*Ny*Nt)], continuous ground truth f0(x,y,t), discrete ground truth df0[Nx*Ny*(Nx*Ny*Nt)]

%% define ground truth function
% continuous GT
f0 = @(x,y,t) gaussmf(x, [sigma_x centering_x]).*gaussmf(y, [sigma_y centering_y]).*gaussmf(t, [sigma_t centering_t]); %continuous ground truth f0(x,t)

%Discrete GT df0
%define number of frames and size of image. Here adapted to gaussian
%parameters.
Frames = centering_t*2+1;
Nt = Frames;
Nx = centering_x*2+1; %size in nb of pixels, x direction
Ny = centering_y*2+1; %size in nb of pixels, y direction
N = Nx*Ny; % size in space
sampling_time = 1/(Nx*Ny);%time to sample 1 pixel

%sample continuous f0
dx = (0:1:Nx-1)'; 
dy = (0:1:Ny-1)'; 
dt(1,1,:) = (0:sampling_time:(N*Frames-1)*sampling_time);

%sampled f0. each element of 3rd dimension is the image in a certain timeframe
%in a certain timeframe x goes from up to down and y from left to right
df0 = f0(repmat(dx,[1,Ny,length(dt)]),repmat(dy',[Nx,1,length(dt)]),repmat(dt,[Nx,Ny,1])); 

%take samples from ground truth
%measurements are considered to be taken first in y direction, then x
%each sample is taken a sampling_time after the previous one
%-------
%|---->|
%|---->|
%|---->|
%-------
measurements = zeros(1,N*Nt);
for i = 1:N*Nt
    measurements(i) = f0(mod(floor(i/Nx),Nx),mod(i,Ny)-1,dt(i));
end

end
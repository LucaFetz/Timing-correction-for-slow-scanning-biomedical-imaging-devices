%% clear workspace
clc
clear all
close all

%% Ground truth simulation
%continuous grounf truth f0
% define size and time ranges
%to do
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
df0 = f0(repmat(dx',1,1),repmat(dt,1,Nx)); %sampled f0. each row is the image in a certain timeframe
df1 = f0(repmat(dx',Nx*Frames,Frames),repmat(dt,1,Nx*Frames)); %sampled f0. In reality only the first Nx columns are needed. This is used to simplify measurement with diag(df1)

%% measurement and usual approximation
measurement = diag(df1); %each sample is taken from its timeframe. In usual approximation it is considered as image at time t=0
figure
plot(repmat(dx,Frames,1), [measurement df1(1,:)' df1(floor(end/2),:)']);
legend(' superposition of Measurements / usual approximation','Ground truth at t=0','Ground truth at time t=1/2')
mean_measurement = mean(reshape(measurement,Nx,Frames),2);
figure
plot(dx, [mean_measurement df1(1,1:Nx)' df1(floor(end/2),1:Nx)']);
legend(' Mean of Measurements / usual approximation','Ground truth at t=0','Ground truth at time t=1/2')

figure 
for i=1:Frames
    plot(dx, [measurement((i-1)*Nx+1:i*Nx) df1(i,(i-1)*Nx+1:i*Nx)']);
    legend(' Measurement / usual approximation','Ground truth')
    title(i)
    pause(0.001)
end

%% forward model
%continuous forward model
k = repmat(0:1:Nx-1,Nt,1);
l = repmat((0:1:Nt-1)',1,Nx);
h = @(x,t) (B_spline(t-l).*B_spline(x-k));
H = zeros(Nx*Nt,Nx*Nt);
for i = 0:1:Nx-1
    for j = 0:1:Nt-1 
       H(i*(Nx-1)+(j+1),:) = reshape(h(i,j+i*sampling_time),1,Nx*Nt);
    end
end
y = diag(measurement);

%% check conditionement of H
%visualize
imagesc(H);
condition = cond(H); %bad condition
%% coefficient optimization for B spline interpolation
%naive inverse approach
C = H\measurement; %inv(H)*measurement;
C = reshape(C,Nt,Nx);
imdisp(C,'C found by Naive approach',1);
%optimization approach

%Makes sense to me but commented because dimensions don't match
% index = 1;
% H2 = LinOpBroadcastMatrix(H,size(measurement),index);  %not really sure how to put H to the good format
% %H2 = LinOpDiag(size(H)); %test. not a good H
% LS=CostL2([],measurement);  % Least-Sqaures data term
% F=LS*H2; %composition of cost and H
% 
% GD=OptiGradDsct(F);
% GD.OutOp=OutputOpti(1,df0,40);
% GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
% GD.maxiter=200;         % max number of iterations
% GD.run(measurement); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

%Dimensions match but not sure it has any meaning
index = 1;
H2 = LinOpBroadcastMatrix(H,size(measurement),index);  %not really sure how to put H to the good format
%H2 = LinOpDiag(size(H)); %test. not a good H
LS=CostL2([],measurement);  % Least-Sqaures data term
F=LS*H2; %composition of cost and H

GD=OptiGradDsct(F);
GD.OutOp=OutputOpti(1,df0(:,1),40); %Why take only part of ground truth df0?? Dimension match but no meaning
GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
GD.maxiter=200;         % max number of iterations
GD.run(measurement); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

C2 = reshape(GD.xopt,Nt,Nx);
%% Display
imdisp(C2,'C found by Gradient Descent',1);

%% reconstruction
reconst = H*reshape(C2,Nx*Nt,1); 
reconst = reshape(reconst,Nt,Nx);
imagesc(reconst)
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
df0 = f0(dx,repmat(dt',Nx,1)); %sampled f0. each column is the image in a certain timeframe
df1 = repmat(df0,Nt,1); %sampled f0. In reality only the first Nx columns are needed. This is used to simplify measurement with diag(df1)
df2 = f0(dx,1:Nt);
%% measurement and usual approximation
measurement = diag(df1); %each sample is taken from its timeframe. In usual approximation it is considered as image at time t=0
figure
plot(repmat(dx,Frames,1), [measurement reshape(df0(:,1:Nt:Nx*(Nt-1)),Nx*Nt,1) reshape(df0(:,(1:Nt:Nx*(Nt-1))+floor(Nt/2)),Nx*Nt,1)]);
legend(' superposition of Measurements / usual approximation','Ground truth at t=0','Ground truth at time t=1/2')
mean_measurement = mean(reshape(measurement,Nx,Frames),2);
figure
plot(dx, [mean_measurement mean(df0(:,1:Nt:Nx*(Nt-1)),2) mean(df0(:,(1:Nt:Nx*(Nt-1))+floor(Nt/2)),2)]);
legend(' Mean of Measurements / usual approximation','Mean of Ground truth at t=integer','Mean of Ground truth at time t=integer+1/2')

figure 
for i=1:Frames
    plot(dx, [measurement((i-1)*Nx+1:i*Nx) df0(:,(i-1)*Nt+1)]); %measurement of image i, ground truth at time i
    legend(' Measurement / usual approximation','Ground truth')
    title(i)
    pause(0.01)
end

%% forward model
%continuous forward model
k = repmat(0:1:Nx-1,Nt,1);
l = repmat((0:1:Nt-1)',1,Nx);
%h rows are fixed in time and columns are fixed in space. represent result
%of B_splines for a given pixel
h = @(x,t) (B_spline(x-k).*B_spline(t-l));%takes scalar input and gives matrix output (h(x-0,t-0), h(x-1,t-0),..h(x-(Nx-1),t-(Nt-1))

H = zeros(Nx*Nt,Nx*Nt); % each row of H is h unfolded.
index = 0;
for j = 0:1:Nt-1 
    for i = 0:1:Nx-1
        %the reshape on h is made to keep the rows of h(x,t) intact and
        %after the reshape we have h(x,t)=[h(x-0,t-0) h(x-1,t-0) h(x-2,t-0)...
        %.. h(x-(Nx-1),t-0) h(x-0,t-1) h(x-1,t-1) ... h(x-(Nx-1),t-(Nt-1))]
        index = index+1;
       H(index,:) = reshape(h(i,j+i*sampling_time)',1,Nx*Nt); %is it better to have the double diag? can be done keeping the columns intact (juste remove')
%rows of H are:
%H=[h(0,0);h(1,sampling_time);h(2,2*sampling_time);...;h((Nx-1),(Nx-1)*sampling_times);
%h(1,1+sampling_time);h(2,1+2*sampling_time);...;h((Nx-1),(Nx-1)*sampling_time+(Nt-1))]
    end
end
y = diag(measurement); %

%% check conditionement of H
%visualize
figure
imagesc(H);
condition = cond(H); %bad condition
%% coefficient optimization for B spline interpolation
%naive inverse approach
C = H\measurement; %inv(H)*measurement;
C = reshape(C,Nx,Nt)';
imdisp(C,'C found by Naive approach',1);
%optimization approach

%Makes sense to me but commented because dimensions don't match
H2 = LinOpMatrix(H);  %not really sure how to put H to the good format
LS=CostL2([],measurement);  % Least-Sqaures data term
F=LS*H2; %composition of cost and H

GD=OptiGradDsct(F);
GD.OutOp=OutputOpti(1,[],40); %prendre df0 avec un pas de 1 et pas dt
GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
GD.maxiter=200;         % max number of iterations
GD.run(measurement); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

%Dimensions match but not sure it has any meaning

% H2 = LinOpBroadcastMatrix(H);  %not really sure how to put H to the good format
% %H2 = LinOpDiag(size(H)); %test. not a good H
% LS=CostL2([],measurement);  % Least-Sqaures data term
% F=LS*H2; %composition of cost and H
% 
lambda = 0*3e-5; %thuner la r�gul
R = CostL2(size(measurement)); %regulation term. should it be a Cost? plut�t size(C)
F2 = F+lambda*R; %how to integrate the regulation term?
%regarder evolcost pour voir nb iterations

% 
% GD=OptiGradDsct(F);
% GD.OutOp=OutputOpti(1,df0(:,1),40); %Why take only part of ground truth df0?? Dimension match but no meaning
% GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
% GD.maxiter=200;         % max number of iterations
% GD.run(measurement); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

C2 = reshape(GD.xopt,Nx,Nt)';
%% Display
imdisp(C2,'C found by Gradient Descent',1);


%% reconstruction
reconst = H*reshape(C2,Nx*Nt,1); 
reconst = reshape(reconst,Nt,Nx);
imagesc(reconst)

%%

figure;
subplot(221); imagesc(reshape(measurement,[11,10]));
axis image;title('Measurements / naive approach');colorbar;
subplot(222); imagesc(df0(:,1:11:end));axis image;title('GT at integer time');colorbar;
subplot(223); imagesc(C2');axis image;
title('Recovered coefficients / Recon (interpolant)');colorbar;
subplot(224); imagesc(df0(:,1:11:end) - C2');axis image;title('GT - Coef');colorbar;

% %% clear workspace and add library
% clc
% clear all
% close all
% path = genpath('GlobalBioIm-release');
% addpath(path);
% 
% %% define parameters for the simulation
% spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
% noise = 0; % 0 or 1
% snr = 10; %signal to noise ratio
% lambda = 0*3e-2; %thune regulation term. 0 for no regulation.

%% take parameters already in workspace
spline_order = param.spline_order; %1 or 3. any other number will result in cubic splines (order 3)
noise = param.noise; % 0 or 1
snr = param.snr; %signal to noise ratio
lambda = param.lambda; %thune regulation term. 0 for no regulation.

%% Ground truth simulation
%continuous ground truth f0
% define centering and sigma of gaussians in space and time
sigma_x = 2;
centering_x = 5;
sigma_t = 2;
centering_t = 5;
% define ground truth function
f0 = @(x,t) gaussmf(x, [sigma_x centering_x]).*gaussmf(t, [sigma_t centering_t]); %continuous ground truth f0(x,t)

%Discrete ground truth df0
%define number of frames and size of image. Here adapted to gaussian
%parameters.
Frames = centering_t*2+1;
Nt = Frames;
Nx = centering_x*2+1; %size in nb of pixels
sampling_time = 1/(Nx);
%sample continuous f0
dx = (0:1:Nx-1)'; 
dt = (0:sampling_time:(Nx*Frames-1)*sampling_time)';

df0 = f0(dx,repmat(dt',Nx,1)); %sampled f0. each column is the image in a certain timeframe
df1 = repmat(df0,Nt,1); %sampled f0. In reality only the first Nx columns are needed. This is used to simplify measurement with diag(df1)
df2 = f0(dx,repmat(1:Nt,Nx,1)); % sampled f0 at integers. why df0(:,1:Nx:end) - df2 !=0??

%% measurement and usual approximation. noise simulation.
measurement = diag(df1); %each sample is taken from its timeframe. In usual approximation it is considered as image at time t=0
noisy_measurement = awgn(measurement,snr,'measured');
if noise == 1
    figure(50), plot(noisy_measurement),hold on, plot(measurement);
    figure(51),subplot(121), surf(reshape(noisy_measurement,Nx,Nt)),subplot(122), surf(reshape(measurement,Nx,Nt));
    measurement = noisy_measurement;
    hold off
end


figure
plot(repmat(dx,Frames,1), [measurement reshape(df0(:,1:Nx:Nt*Nx),Nx*Nt,1) reshape(df0(:,((1:Nx:Nx*Nt)+floor(Nt/2))),Nx*Nt,1)]);
legend(' superposition of Measurements / usual approximation','Ground truth at t=integer','Ground truth at time t=integer+1/2')
mean_measurement = mean(reshape(measurement,Nx,Frames),2);
figure
plot(dx, [mean_measurement mean(df0(:,1:Nx:Nt*Nx),2) mean(df0(:,(1:Nx:Nt*Nx)+floor(Nt/2)),2)]);
legend(' Mean of Measurements / usual approximation','Mean of Ground truth at t=integer','Mean of Ground truth at time t=integer+1/2')

%slice animation
% figure 
% for i=1:Frames
%     plot(dx, [measurement((i-1)*Nx+1:i*Nx) df0(:,i*Nx)]); %measurement of image i, ground truth at time i
%     legend(' Measurement / usual approximation','Ground truth')
%     title(i)
%     pause(0.01)
% end

%% forward model
%continuous forward model
k = repmat(0:1:Nx-1,Nt,1);
l = repmat((0:1:Nt-1)',1,Nx);
%h rows are fixed in time and columns are fixed in space. represent result
%of (cubic_)B_splines for a given pixel
if spline_order == 1
    h = @(x,t) (B_spline(x-k).*B_spline(t-l));%takes scalar input and gives matrix output (h(x-0,t-0), h(x-1,t-0),..h(x-(Nx-1),t-(Nt-1))
else
    h = @(x,t) (cubic_B_spline(x-k).*cubic_B_spline(t-l));%takes scalar input and gives matrix output (h(x-0,t-0), h(x-1,t-0),..h(x-(Nx-1),t-(Nt-1))
end

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


%% check conditionement of H
%visualize
figure
imagesc(H);
title('Forward model H')
condition = cond(H); %bad condition on peut am�liorer avec H+lambda(eye)

%% coefficient optimization for B spline interpolation
%naive inverse approach
C = H\measurement; %inv(H)*measurement;
C = reshape(C,Nx,Nt)';
figure, imagesc(C);
title('C found by Naive approach');

%optimization approach
H2 = LinOpMatrix(H);  %put H to LinOpMatrix format
LS=CostL2([],measurement);  % Least-Sqaures data term
F=LS*H2; %composition of cost and H
R = CostL2(size(C(:))); %regulation term
F2 = F+lambda*R; %how to integrate the regulation term?
%regarder evolcost pour voir nb iterations

GD=OptiGradDsct(F2); %GradientDescent
GD.OutOp=OutputOpti(1,[],40); %prendre df0 avec un pas de 1 et pas dt
GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
GD.maxiter=500;         % max number of iterations
GD.run(measurement); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

C2 = reshape(GD.xopt,Nx,Nt)';

%% Direct inversion (show numerical instability that occurs rapidly)
lambda2 = 1e0;
condition = cond(H  + lambda2*eye(size(H)));
C = (H  + lambda2*eye(size(H)))\measurement;
figure;clf
imagesc(reshape(C,Nx,Nt));axis image;colorbar;
title('C found by direct inversion')
F*C;
%% Display
figure, imagesc(C2);
title('C found by Gradient Descent');

%% reconstruction of measurements from model
reconst = H2*GD.xopt;%H*C2 properly shaped gives back measurement
reconst = reshape(reconst,Nx,Nt);

%compares measurement and H*C
figure;
subplot(221);imagesc(reconst);axis image;colorbar;
title('reconstruction of measurements');
subplot(222), imagesc(reshape(measurement,Nx,Nt));axis image;colorbar;
title('measurements')
subplot(223);imagesc(reconst - reshape(measurement,Nx,Nt));axis image;colorbar;
title('reconstruction - measurements')

%% comparison of measurements, GT and recovered coefficients

figure;
subplot(221); imagesc(reshape(measurement,[Nx,Nt]));
axis image;title('Measurements / naive approach');colorbar;caxis([0,max(df0(:))]);

subplot(222); imagesc(df0(:,1:Nx:end));axis image;
title('GT at integer time');colorbar;caxis([0,max(df0(:))]);

subplot(223); imagesc(C2');axis image;
title('Recovered coefficients / Recon (interpolant)');colorbar;caxis([0,max(df0(:))]);

subplot(224); imagesc(df0(:,1:Nx:end) - C2');axis image;title('GT - Coef');
colorbar;
linkaxes(findall(gcf,'Type','Axes'),'xy');

%% interpolating function f

f = @(x,t) sum(sum(C2.*h(x,t)));
df = zeros(Nx,Nx*Nt);
%samling f in df with double loop. maybe try to vectorize
index = 0;
for j = (0:sampling_time:(Nx*Frames-1)*sampling_time)
    index = index +1;
    for i = 0:1:Nx-1
        df(i+1,index) = f(i,j);
    end
end
%df = f(dx,repmat(dt',Nx,1)); %sampled f0. each column is the image in a certain timeframe
%comparison between sampled f and sampled GT in 2D
figure, subplot(221),imagesc(df);
title('Sampled f(x,t) / reconstructed image')
subplot(222),imagesc(df0);
title('Sampled GT')
subplot(223),imagesc(df-df0);

%comparison between sampled f, GT (at integer times) and measures in 3D
figure, subplot(231),surf(df0(:,1:Nx:Nx*Nt));zlim([0 max(max(df(:,1:Nx:Nx*Nt)))]);
title('Sampled GT')
subplot(232),surf(df(:,1:Nx:Nx*Nt));zlim([0 max(max(df(:,1:Nx:Nx*Nt)))]);
title('Corrected measurements')
subplot(233),surf(reshape(measurement,Nx,Nt));zlim([0 max(max(df(:,1:Nx:Nx*Nt)))]);
title('Measurements/Naive approach')
subplot(235),surf(df0(:,1:Nx:Nx*Nt) - df(:,1:Nx:Nx*Nt));
title('GT - Corrected measurements')
subplot(236),surf(df0(:,1:Nx:Nx*Nt) - reshape(measurement,Nx,Nt));
title('GT - Measurements/Naive approach')

%comparison between continuous f and GT
figure, subplot(211),fsurf(f0,[0 Nx-1 0 Nt-1]);
title('continuous GT')
subplot(212),fsurf(f,[0 Nx-1 0 Nt-1]);
title('continuous f(x,t)')


%% Alternative way to get the samples

% ind = [1+mod(0:Nx*Nt-1,Nx)',(1:Nx*Nt)'];
% tmp_samples = df0;
% for kk = 1:size(ind,1)
% tmp_samples(ind(kk,1),ind(kk,2)) = inf;
% end
% figure,imagesc(tmp_samples),colorbar;
% xticks(1:20:Nx*Nt);
% xticklabels(0:20:Nx*Nt-1), yticklabels(0:Nx-1);
% xlabel('Time [dt]'),ylabel('Space [px]');
% title('Location of samples on ground truth')
% 
% tmp_samples = df0;
% tmp_samples(:,1:Nx:Nx*Nt) = inf;
% figure,imagesc(tmp_samples),colorbar;
% xticks(1:20:Nx*Nt);
% xticklabels(0:20:Nx*Nt-1), yticklabels(0:Nx-1);
% xlabel('Time [dt]'),ylabel('Space [px]');
% title('Samples on ground truth at time t')
% 
% 
% tmp_samples = zeros(size(df0));
% for i = 1:size(df0,1)
%         tmp_samples(i,1:(end-i+1)) = df0(i,i:end);
% end
% 
% tmp_samples(:,1:Nx:Nx*Nt) = inf;
% figure,imagesc(tmp_samples),colorbar;
% xticks(1:20:Nx*Nt);
% xticklabels(0:20:Nx*Nt-1), yticklabels(0:Nx-1);
% xlabel('Time [dt]'),ylabel('Space [px]');
% title('Distorsion on ground truth due to approximation')
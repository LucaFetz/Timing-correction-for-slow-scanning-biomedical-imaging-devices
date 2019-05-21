function [measurements, f0, GT, Nx, Ny, Nt, samples_coordinates_x, samples_coordinates_y] = create_2d_ground_truth(samples_coordinates_x, samples_coordinates_y, sigma_x,sigma_y,sigma_t,centering_x,centering_y, centering_t, noise, noise_snr, config, speed)
%function that creates 2D+time ground truth(GT) as a 2D gaussian in space and a
%gaussian in time. Takes samples from GT and adds noise if needed.
%input: vector defining position of samples. Can be 'classic' or
%'backnforth' to use predefined scanning behaviours
%, sigma and centering to define the gaussians
%output:measurements[(Nx*Ny*Nt)], continuous ground truth f0(x,y,t), discrete ground truth GT at integer frames[(Nx*Ny*Nt)]

%% define ground truth function
% continuous GT
switch config
    case 1 %2D temporal gaussian
        f0 = @(x,y,t) gaussmf(x, [sigma_x centering_x]).*gaussmf(y, [sigma_y centering_y]).*gaussmf(t, [sigma_t centering_t]); %continuous ground truth f0(x,t)
    case 2 %multiple 2D temporal gaussian
        f0 = @(x,y,t) gaussmf(x, [sigma_x/2 centering_x/2]).*gaussmf(y, [sigma_y/2 centering_y/2]).*gaussmf(t, [sigma_t/2 centering_t/2]) + gaussmf(x, [sigma_x centering_x]).*gaussmf(y, [sigma_y*1.5 centering_y*1.5]).*gaussmf(t, [sigma_t/2 centering_t*1.5]) + gaussmf(x, [sigma_x centering_x]).*gaussmf(y, [sigma_y centering_y]).*gaussmf(t, [sigma_t centering_t]); %continuous ground truth f0(x,t)
    case 3 %moving gaussian
        r = min(centering_x,centering_y)*0.66; 
        theta = speed/(r);
        f0 = @(x,y,t) gaussmf(x, [sigma_x/2 r * cos(theta*t) + centering_x]).*gaussmf(y, [sigma_y/2 r * sin(theta*t) + centering_y]);
    case 4 %"beating" gaussian
        a1 = 0.5 * sigma_x;
        sigma_x_beating = @(t) a1*sin(speed*t) + sigma_x;
        a2 = 0.5 * sigma_y;
        sigma_y_beating = @(t) a2*sin(speed*t) + sigma_y;
        f0 = @(x,y,t) gaussmf(x, [sigma_x_beating(t) centering_x]).*gaussmf(y, [sigma_y_beating(t) centering_y]);
end
%Discrete GT df0
%define number of frames and size of image. Here adapted to gaussian
%parameters.
Frames = centering_t*2+1;
Nt = Frames;
Nx = centering_x*2+1; %size in nb of pixels, x direction
Ny = centering_y*2+1; %size in nb of pixels, y direction
N = Nx*Ny; % size in space
sampling_time = 1/(Nx*Ny);%time to sample 1 pixel

%% sample continuous f0
dx = (0:1:Nx-1)'; 
dy = (0:1:Ny-1)'; 
dt(1,1,:) = (0:sampling_time:(N*Frames-1)*sampling_time);
 
%take grount truth frames as if acquisition was instantaneous
GT = zeros(1,N*Nt);
for i = 1:N*Nt
    GT(i) = f0(mod(floor((i-1)/Ny),Nx),mod(i-1,Ny),floor(dt(i)));
end
GT = permute(reshape(GT,[Ny,Nx,Nt]),[2 1 3]);

%% --- take samples from ground truth at coordinates specified by user ---
measurements = zeros(1,N*Nt);
%classic scanning behaviour
if(samples_coordinates_x == "classic")
    %in classic measurements are considered to be taken first in y direction, then x
    %each sample is taken a sampling_time after the previous one
    %-------
    %|---->|
    %|---->|
    %|---->|
    %-------
    samples_coordinates_x = mod(floor((0:1:N*Nt-1)/Ny),Nx);
else

    %backnforth scanning behaviour
    if(samples_coordinates_x == "backnforth")
        %in backnforth measurements are considered to be taken first in y
        %direction, then x. At next frame, they are taken in opposite
        %direction.(from bottom right to top left)
        %each sample is taken a sampling_time after the previous one
        %-------        %-------
        %|---->|        %|<---||
        %|---->|    ->  %|<---||    
        %|---->|        %|<---||
        %-------        %-------
        samples_coordinates_x = repmat([repelem(0:Nx-1,Ny), repelem(Nx-1:-1:0,Ny)],1,floor(Nt/2));
        if(floor(Nt/2) ~= Nt/2)
            samples_coordinates_x = [samples_coordinates_x, repelem(0:Nx-1,Ny)];
        end
    end
end
if(samples_coordinates_y == "classic")
    %in classic measurements are considered to be taken first in y direction, then x
    %each sample is taken a sampling_time after the previous one
    %-------
    %|---->|
    %|---->|
    %|---->|
    %-------
    samples_coordinates_y = mod((0:1:N*Nt-1),Ny);
else

    if(samples_coordinates_y == "backnforth")
        %in backnforth measurements are considered to be taken first in y
        %direction, then x. At next frame, they are taken in opposite
        %direction.(from bottom right to top left)
        %each sample is taken a sampling_time after the previous one
        %-------        %-------
        %|---->|        %|<---||
        %|---->|    ->  %|<---||    
        %|---->|        %|<---||
        %-------        %-------
        samples_coordinates_y = repmat([repmat(0:1:Ny-1,[1,Nx]), repmat(Ny-1:-1:0,[1,Nx])],1,floor(Nt/2));

        if(floor(Nt/2) ~=Nt/2)
            samples_coordinates_y = [samples_coordinates_y, repmat(0:1:Ny-1,[1,Nx])];
        end
    end
end

  
% take samples at given coordinates
for i = 1:N*Nt
    measurements(i) = f0(samples_coordinates_x(i),samples_coordinates_y(i),dt(i));
end

%old way to take measurements
% for i = 1:N*Nt
%     measurements(i) = f0(mod(floor((i-1)/Ny),Nx),mod(i-1,Ny),dt(i));
% end

%add white gaussian noise if needed
if noise == 1
    measurements = awgn(measurements,noise_snr,'measured');
end

end
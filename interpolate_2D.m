function [f,reconstructed_frames] = interpolate_2D(h, C)
%create an interpolating function f(x,y,t) from B splines and already
%optimized coefficients.
%input: continuous forward model h, coefficients C
%output: interpolating function f, sampled f: df

%interpolating function
f = @(x,y,t) sum(sum(sum(C.*h(x,y,t))));

%get back sizes of simulation
[Nx, Ny, Nt] = size(h(0,0,0));
sampling_time = 1/(Nx*Ny);
dt(1,1,:) = (0:sampling_time:(Nx*Ny*Nt-1)*sampling_time);
%df = zeros(Nx,Ny,Nx*Ny*Nt);

%sampling reconstructed frames from f with triple loop. maybe try to vectorize
reconstructed_frames = zeros(1,Nx*Ny*Nt);
for i = 1:Nx*Ny*Nt
    reconstructed_frames(i) = f(mod(floor((i-1)/Ny),Nx),mod(i-1,Ny),floor(dt(i))); %sample like measurements but at integer times
end
%reshape the frames
reconstructed_frames = permute(reshape(reconstructed_frames,[Ny,Nx,Nt]),[2 1 3]);
end
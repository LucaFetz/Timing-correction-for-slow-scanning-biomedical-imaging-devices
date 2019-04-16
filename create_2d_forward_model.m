function [H, h] = create_2d_forward_model(Nx,Ny,Nt,spline_order,samples_coordinates_x, samples_coordinates_y)
%Function that returns forward model for sampling 2d images
%input: Nx,Ny,Nt size of images. spline_order, 1 or 3
%Output: forward model H

%% continuous forward model
sampling_time = 1/(Nx*Ny);
k = repmat((0:1:Nx-1)',[1,Ny,Nt]);
l = repmat((0:1:Ny-1),[Nx,1,Nt]);
dt(1,1,:) = (0:1:Nt-1);
m = repmat(dt,[Nx,Ny,1]);
%h 3rd dimension is fixed in time
%h rows are fixed in x and columns are fixed in y. represent result
%of (cubic_)B_splines for a given pixel
if spline_order == 1
    h = @(x,y,t) (B_spline(x-k).*B_spline(y-l).*B_spline(t-m));%takes scalar input and gives matrix output (h(x-0,t-0), h(x-1,t-0),..h(x-(Nx-1),t-(Nt-1))
else
    h = @(x,y,t) (cubic_B_spline(x-k).*cubic_B_spline(y-l).*cubic_B_spline(t-m));%takes scalar input and gives matrix output (h(x-0,t-0), h(x-1,t-0),..h(x-(Nx-1),t-(Nt-1))
end

%% discrete forward model H
H = zeros(Nx*Ny*Nt,Nx*Ny*Nt); % each row of H is h unfolded.
index = 0;

dt2(1,1,:) = (0:sampling_time:(Nx*Ny*Nt-1)*sampling_time);
samples_coordinates_f = repelem(0:Nt-1,1,Nx*Ny);
for k = 1:Nx*Ny*Nt
    h1 = h(samples_coordinates_x(k),samples_coordinates_y(k),dt2(k));
    for l = 1:Nx*Ny*Nt
        h2(l) = h1(samples_coordinates_x(l)+1,samples_coordinates_y(l)+1,samples_coordinates_f(l)+1);
    end
    H(k,:) = h2; 
end

% for k = 1:Nx*Ny*Nt
%     H(k,:) = reshape(permute(h(samples_coordinates_x(k),samples_coordinates_y(k),dt2(k)),[2 1 3]),[1,Nx*Ny*Nt]); 
% end

% old version compatible only with classic scanning
% for k = 0:1:Nt-1
%     for i = 0:1:Nx-1
%         for j = 0:1:Ny-1 
%             %the reshape on h is made to keep the rows of h(x,t) intact, then the columns and
%             %then the 3rd dimension.after the reshape we have h(x,t)=[h(x-0,y-0,t-0) h(x-0,y-1,t-0) h(x-0,y-2,t-0)...
%             %.. h(x-0,y-(Ny-1),t-0) h(x-1,y-0,t-0) h(x-1,y-1,t-0) ...
%             %h(x-(Nx-1),y-(Ny-1),t-0] h(x-0,y-0,t-1]... h(x-(Nx-1),y-(Ny-1),t-(Nt-1)] 
%             index = index+1;
%             H(index,:) = reshape(permute(h(i,j,(index-1)*sampling_time),[2 1 3]),[1,Nx*Ny*Nt]); 
%             %rows of H are:
%             %H=[h(0,0);h(1,sampling_time);h(2,2*sampling_time);...;h((Nx-1),(Nx-1)*sampling_times);
%             %h(1,1+sampling_time);h(2,1+2*sampling_time);...;h((Nx-1),(Nx-1)*sampling_time+(Nt-1))]
%         end
%     end
% end
%since H is full of zeros, store it correctly
H = sparse(H);
end

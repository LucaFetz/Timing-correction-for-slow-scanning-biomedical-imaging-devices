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

dt2(1,1,:) = (0:sampling_time:(Nx*Ny*Nt-1)*sampling_time);
samples_coordinates_f = repelem(0:Nt-1,1,Nx*Ny);
h2 = zeros(1,Nx*Ny*Nt);
for k = 1:Nx*Ny*Nt
    h1 = h(samples_coordinates_x(k),samples_coordinates_y(k),dt2(k));
    for l = 1:Nx*Ny*Nt
        h2(l) = h1(samples_coordinates_x(l)+1,samples_coordinates_y(l)+1,samples_coordinates_f(l)+1);
    end
    H(k,:) = h2; 
end
%since H is full of zeros, store it correctly
H = sparse(H);

%% work in progress, optimize H creation filling only nonzero members of sparse matrix
% -- continuous h --
% sampling_time = 1/(Nx*Ny);
% dt2(1,1,:) = (0:sampling_time:(Nx*Ny*Nt-1)*sampling_time);
% samples_coordinates_f = repelem(0:Nt-1,1,Nx*Ny);
% %h gives B_splines result for a given pixel x,y,t for given coefficient
% % coordinates k,l,m
% if spline_order == 1
%     h2 = @(x,y,t,k,l,m) (B_spline(x-k).*B_spline(y-l).*B_spline(t-m));%takes scalar input and gives scalar output
% else
%     h2 = @(x,y,t,k,l,m) (cubic_B_spline(x-k).*cubic_B_spline(y-l).*cubic_B_spline(t-m));%takes scalar input and gives scalar output
% end
% 
% % -- discrete H --
% 
% %Hopt = sparse(Nx*Ny*Nt,Nx*Ny*Nt);
% %Hopt2 = sparse(Nx*Ny*Nt,Nx*Ny*Nt);
% %Hopt3 = sparse(Nx*Ny*Nt,Nx*Ny*Nt);
% tmp = 0;
% i = ones(Nx*Ny*Nt*3*3*4,1);
% j = ones(Nx*Ny*Nt*3*3*4,1);
% v = zeros(Nx*Ny*Nt*3*3*4,1);
% for k = 1:Nx*Ny*Nt
%     %h1 = zeros(Nx,Ny,Nt);
%     
%     %triple loop to fill only the non zero coeffs
%     for k2 = 0:2
%         for l2 = 0:2
%             for m2 = 1:4
%                 %pixel order in H follows the scanning path
%                 idx1 = samples_coordinates_x(k)+k2;
%                 idx2 = samples_coordinates_y(k)+l2;
%                 idx3 = floor(dt2(k))-1+m2;
%                 
%                 %avoid borders
%                 if(idx1 > 0 && idx2 > 0 && idx3 > 0 && idx1 <= Nx && idx2 <= Ny && idx3 <= Nt)
%                     %h1(idx1,idx2,idx3) = h2(samples_coordinates_x(k),samples_coordinates_y(k),dt2(k),samples_coordinates_x(k)-1+k2,samples_coordinates_y(k)-1+l2,floor(dt2(k))-2+m2); %remplir h1 avec boucle
%                     
%                     %Hopt2(k,(idx1-1)*Ny + idx2 + (idx3-1)*Nx*Ny) = h2(samples_coordinates_x(k),samples_coordinates_y(k),dt2(k),samples_coordinates_x(k)-1+k2,samples_coordinates_y(k)-1+l2,floor(dt2(k))-2+m2);%remplir Hopt2 directement avec h
%                     
%                     tmp = tmp+1;
%                     i(tmp) = k;
%                     j(tmp) = (idx1-1)*Ny + idx2 + (idx3-1)*Nx*Ny;
%                     v(tmp) = h2(samples_coordinates_x(k),samples_coordinates_y(k),dt2(k),samples_coordinates_x(k)-1+k2,samples_coordinates_y(k)-1+l2,floor(dt2(k))-2+m2);
%                     
%                 end
%             end
%         end
%     end
%     %each row of H has to follow scannig path order
%     for l = 1:Nx*Ny*Nt
%         %Hopt(k,l) = h1(samples_coordinates_x(l)+1,samples_coordinates_y(l)+1,samples_coordinates_f(l)+1);         
%     end
%    
% end
% %each row of Hopt/Hopt3 has the coeffs for 1 pixel of the image, in
% %scanning path order. the rows are also in scanning path order.
% 
% Hopt4 = sparse(i,j,v,Nx*Ny*Nt,Nx*Ny*Nt);
% Hopt5 = sparse(size(Hopt4));
% 
% for l = 1:Nx*Ny*Nt
%     %Hopt3(:,l) = Hopt2(:,samples_coordinates_x(l)*Ny+samples_coordinates_y(l)+1+samples_coordinates_f(l)*Ny*Nx);
% 
%     Hopt5(:,l) = Hopt4(:,samples_coordinates_x(l)*Ny+samples_coordinates_y(l)+1+samples_coordinates_f(l)*Ny*Nx);
% 
% end
breakpoint = 0;

%% old version compatible only with classic scanning
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
%H = sparse(H);
end

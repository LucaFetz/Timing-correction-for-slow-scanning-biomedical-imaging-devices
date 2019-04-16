function result = Etape_2D(param)

%% Create 2D ground truth
%define GT parameters
param.GT.sigma_x = 2;
param.GT.centering_x = 3;
param.GT.sigma_y = 2;
param.GT.centering_y = 4;
param.GT.sigma_t = 2;
param.GT.centering_t = 5;
%create GT
[measurements, f0, GT, Nx, Ny, Nt] = create_2d_ground_truth(param.GT.sigma_x,param.GT.sigma_y,param.GT.sigma_t,param.GT.centering_x,param.GT.centering_y, param.GT.centering_t, param.noise, param.noise_snr);

%% plot GT with animation in time
if(param.plot_flag)
    figure;
    for i = 1:Nt
        imagesc(GT(:,:,i));colorbar;caxis([0 1]);
        title('GT')
        pause(0.2);
    end
end
%% plot measurements with animation in time

    %since the reshape is column wise we reshape as [Ny,Nx,Nt] and then permute
    measurements_matrix = permute(reshape(measurements,[Ny,Nx,Nt]),[2 1 3]);
if(param.plot_flag)
    figure;
    for i = 1:Nt
        imagesc(measurements_matrix(:,:,i));colorbar;caxis([0 1]);
        title('measurements/naive approach')
        pause(0.2);
    end
end
%% Create forward model
[H, h] = create_2d_forward_model(Nx,Ny,Nt,param.spline_order);

%% plot H
if(param.plot_flag)
    figure, imagesc(H);
end
%% Find C with inverse problem
C = optimize_c_2D(H, param.lambda, measurements',param.opti_type,param.regul_type);
C = permute(reshape(C',[Ny,Nx,Nt]),[2 1 3]);

%% plot C
if(param.plot_flag)
    figure
    for i = 1:Nt
        imagesc(C(:,:,i));colorbar;caxis([0 1]);
        title('Coefs');
        pause(0.2);
    end
end
%% Reconstruct good frames
[f,result.reconstructed_frames] = interpolate_2D(h, C);
%% visualize reconstructed frames
if(param.plot_flag)
    figure;
    for i = 1:Nt
        imagesc(result.reconstructed_frames(:,:,i));colorbar;caxis([0 1]);
        title('reconstructed frames');
        pause(0.2);
    end
end
%% assess reconstruction quality and visualize results


    result.snr_measurements = snr(GT,GT-measurements_matrix);
    result.snr_reconstruction = snr(GT,GT-result.reconstructed_frames);
if(param.plot_flag)
    figure('Units','normalized','Position',[0 0 1 1]);

    for i = 1:Nt

        %str = sprintf('Time = %d', i)
        %suptitle(str)

        subplot(231),imagesc(GT(:,:,i));colorbar;caxis([0 1]);
        title('GT')

        subplot(232),imagesc(measurements_matrix(:,:,i));colorbar;caxis([0 1]);
        title('measurements/naive approach')

        subplot(233),imagesc(result.reconstructed_frames(:,:,i));colorbar;caxis([0 1]);
        title('reconstructed frames');

        subplot(235),imagesc(GT(:,:,i) - measurements_matrix(:,:,i));colorbar;caxis([0 max(max(max(GT(:,:,:) - measurements_matrix(:,:,:))))]);
        title('GT - measurements');

        subplot(236),imagesc(GT(:,:,i) - result.reconstructed_frames(:,:,i));colorbar;caxis([0 max(max(max(GT(:,:,:) - result.reconstructed_frames(:,:,:))))]);
        title('GT - reconstructed frames');

        pause(0.2);
    end

    for i = 1:Nt

        %str = sprintf('Time = %d', i)
        %suptitle(str)

        subplot(231),surf(GT(:,:,i));colorbar;caxis([0 1]);zlim([0 1]);
        title('GT')

        subplot(232),surf(measurements_matrix(:,:,i));colorbar;caxis([0 1]);zlim([0 1]);
        title('measurements/naive approach')

        subplot(233),surf(result.reconstructed_frames(:,:,i));colorbar;caxis([0 1]);zlim([0 1]);
        title('reconstructed frames');

        subplot(235),surf(GT(:,:,i) - measurements_matrix(:,:,i));colorbar;caxis([min(min(min(GT(:,:,:) - measurements_matrix(:,:,:)))) max(max(max(GT(:,:,:) - measurements_matrix(:,:,:))))]);zlim([min(min(min(GT(:,:,:) - measurements_matrix(:,:,:)))) max(max(max(GT(:,:,:) - measurements_matrix(:,:,:))))]);
        title('GT - measurements');

        subplot(236),surf(GT(:,:,i) - result.reconstructed_frames(:,:,i));colorbar;caxis([min(min(min(GT(:,:,:) - result.reconstructed_frames(:,:,:)))) max(max(max(GT(:,:,:) - result.reconstructed_frames(:,:,:))))]);zlim([min(min(min(GT(:,:,:) - result.reconstructed_frames(:,:,:)))) max(max(max(GT(:,:,:) - result.reconstructed_frames(:,:,:))))]);
        title('GT - reconstructed frames');

        pause(0.2);
    end
end

end
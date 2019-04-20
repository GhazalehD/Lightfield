function depth_reg_image = confi_depth(LF)

addpath(genpath('required'));



% INTERNAL PARAMETERS

%%% LF sizes                        --------------
UV_radius           = 3;
UV_diameter         = (2*UV_radius+1);
UV_size             = UV_diameter^2;

%%% Shearing                        --------------
depth_resolution        = 100;
alpha_min               = 0.2;
alpha_max               = 2;

%%% Analysis                        --------------
% defocus analysis radius
defocus_radius          = 3;
% correspondence analysis radius
corresp_radius          = UV_radius;

%%% Regularize                      --------------
lambda_data             = 1;
lambda_flat             = 2;
lambda_smooth           = 1;
lambda_shading          = 0;
confidence_scale        = 0.6; 
iter_max                = 20;
convergence_ratio       = 0.00001;
err_epsilon             = 0.001;
smooth_kernel           = [0 -1 0;-1 4 -1;0 -1 0];
flat_kernel             = [0 0 0;-1 1 0;0 0 0];

%%% Reprojection
pixel_length            = 0.001;

%%% Shading
tau_r                   = 0.0005;
normals_num_clusters    = 200;
texture_num_clusters    = 200;
k_means_iter_max        = 500;
k_means_spatial_weight  = 0.0001;
max_connections         = 4;    
w_local_shading         = 1;
w_local_reflectance     = 1;
w_nonlocal_shading      = 0.05;
w_nonlocal_reflectance  = 0.01;

    LF = LF(3:9,3:9,:,:,1:3);
    % LOAD CAMERA DATA
    x_size              = size(LF,4);
    y_size              = size(LF,3);
    % JPEG (RAW IMAGE)
    LF_x_size           = size(LF,4)*UV_diameter;
    LF_y_size           = size(LF,3)*UV_diameter;                                                               
    % read out file



% GATHER PARAMTERS
    LF_parameters       = struct('LF_x_size',LF_x_size,...
        'LF_y_size',LF_y_size,...
        'x_size',x_size,...
        'y_size',y_size,...
        'UV_radius',UV_radius,...
        'UV_diameter',UV_diameter,...
        'UV_size',UV_size,...
        'depth_resolution',depth_resolution,...
        'alpha_min',alpha_min,...
        'alpha_max',alpha_max,...
        'defocus_radius',defocus_radius,...
        'corresp_radius',corresp_radius,...
        'lambda_data',lambda_data,...
        'lambda_flat',lambda_flat,...
        'lambda_smooth',lambda_smooth,...
        'lambda_shading',lambda_shading,...
        'confidence_scale',confidence_scale,...
        'iter_max',iter_max,...
        'convergence_ratio',convergence_ratio,...
        'err_epsilon',err_epsilon,...
        'smooth_kernel',smooth_kernel,...
        'flat_kernel',flat_kernel,...
        'pixel_length',pixel_length,...
        'tau_r',tau_r,...
        'normals_num_clusters',normals_num_clusters,...
        'texture_num_clusters',texture_num_clusters,...
        'k_means_iter_max',k_means_iter_max,...
        'k_means_spatial_weight',k_means_spatial_weight,...
        'max_connections',max_connections,...
        'w_local_shading',w_local_shading,...
        'w_local_reflectance',w_local_reflectance,...
        'w_nonlocal_shading',w_nonlocal_shading,...
        'w_nonlocal_reflectance',w_nonlocal_reflectance...
        ) ;

% GATHER NECESSARY DATA
fprintf('I. Remapping LF JPEG to our standard                  *******\n');
tic;
% RAW to Remap                    
       
LF_Remap            = RAW2REMAP1...
                        (LF,LF_parameters)    ;
% Remape to Pinhole (Center View)
IM_Pinhole          = squeeze(LF(4,4,:,:,:));




[defocus_response, corresp_response] = compute_LFdepth(LF_Remap, IM_Pinhole, LF_parameters);

%%

defocus_confi = ALM_CONFIDENCE(defocus_response,IM_Pinhole,LF_parameters,0);
corresp_confi = ALM_CONFIDENCE(corresp_response,IM_Pinhole,LF_parameters,0);

[defocus_confi,corresp_confi] = NORMALIZE_CONFIDENCE(defocus_confi,corresp_confi);     

combined_response = COMBINE_RESPONSES(defocus_response, corresp_response, defocus_confi, corresp_confi);

combined_depth = DEPTH_ESTIMATION(combined_response,0);


combined_confi = ALM_CONFIDENCE(combined_response,IM_Pinhole,LF_parameters,0);
combined_confi = (combined_confi - min(combined_confi(:)))/(max(combined_confi(:)) - min(combined_confi(:)));



[depth_reg_image] = depth_regularization(IM_Pinhole, combined_depth, combined_confi, LF_parameters);

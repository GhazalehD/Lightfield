% Loading Image
clear 
close all
addpath(genpath('matchSIFT'))
lighfield= load('pen.mat');
lighfield =im2double(lighfield.LF(:,:,:,:,1:3));

H =[0.0002         0         0         0   -0.0013
         0    0.0002         0         0   -0.0013
         0         0    0.0004         0   -0.0720
         0         0         0    0.0004   -0.0720
         0         0         0         0    1.0000];
f =     0.1999;

%% Extracting SIFT Features For the Centeral Sub-aperture

edge_thresh = 2; % you will have more SIFT features if you use a higher threshold.
IM_center = rgb2gray(squeeze(lighfield(6,6,:,:,:)));

[SIFT_Loc,SIFT_des] = vl_sift(single(IM_center), 'edgethresh', edge_thresh) ; % The outputs will be SIFT Locations and SIFT decriptors 
SIFT_Loc = SIFT_Loc(1:2,:)'; % The First and second rows store the X, Y locations of each SIFT feature respectively. see vl_sift for more details.

figure(1);imshow(squeeze(lighfield(6,6,:,:,:)))
hold on
plot(SIFT_Loc(:,1),SIFT_Loc(:,2),'*g')

%% SIFT Matching between SIFT From the Centeral Sub-aperture and all other Sub-apertures

tic
N_SIT = size(SIFT_Loc,1);
match = zeros(N_SIT,25*4);
cnt = 1;
% Note that we use only 25 sub-apertures to find the SIFT matching, outer
% sub-aperture are almost useless due to Vignetting problem.
% and also using lesser sub-aprture will yeid more computational speed.
for i = 4:8
    for j = 4:8
        
        gray = rgb2gray(squeeze(lighfield(j,i,:,:,:)));
        % calculating SIFT for each Sub-apertures
        [SIFTloc_j,SIFTdes_j] = vl_sift(single(gray), 'edgethresh', edge_thresh) ;
        
        % finding the indices of corrspoding SIFT features
        [matchi, matchj] = matchSIFTdesImagesBidirectional(SIFT_des, SIFTdes_j);
        
        % Storing Location of SIFT corrspondences with their sub-aperture indices  
        % each row is corresponding to a specific SIFT feature
        match(matchi,cnt:cnt+3)= [repmat([i j],length(matchi),1) SIFTloc_j(1:2,matchj)' ];
        cnt = cnt+4;
    end
end
toc

%% Calculating the 3D Location of each SIFT Features

P3D = zeros(N_SIT,3);
for counter = 1: N_SIT
    
    corres = reshape(match(counter,:),4,25)';
    Nz = find(corres(:,1)~=0); % it might be in some sub-apetures we cannot find correspondencing features. 
    
    % Note that we calculate the 3D location only if the number of correspondences feature are more than 10 
    % less than this number will result in inaccurate 3D Location.
    if length(Nz)>10      
        rays  = H*[corres(Nz,:) ones(length(Nz),1)]'; % transforming pixel coordinates to ray coordinates using intrinsic matrix H
        s = Point_3Dlocation_estimate( rays );
        fz1 =  -s(1)+1;
        P3D(counter,:) = [s(2)/fz1 s(3)/fz1 f/fz1];
    end
end
%%
K = [SIFT_Loc P3D(:,3)];
lighfield= load('pen.mat');
lighfield =im2double(lighfield.LF(:,:,:,:,1:3));

H =[0.0002         0         0         0   -0.0013
         0    0.0002         0         0   -0.0013
         0         0    0.0004         0   -0.0720
         0         0         0    0.0004   -0.0720
         0         0         0         0    1.0000];
f =     0.1999;

% Extracting SIFT Features For the Centeral Sub-aperture

edge_thresh = 2; % you will have more SIFT features if you use a higher threshold.
IM_center = rgb2gray(squeeze(lighfield(6,6,:,:,:)));

[SIFT_Loc,SIFT_des] = vl_sift(single(IM_center), 'edgethresh', edge_thresh) ; % The outputs will be SIFT Locations and SIFT decriptors 
SIFT_Loc = SIFT_Loc(1:2,:)'; % The First and second rows store the X, Y locations of each SIFT feature respectively. see vl_sift for more details.

% SIFT Matching between SIFT From the Centeral Sub-aperture and all other Sub-apertures

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
% Calculating the 3D Location of each SIFT Features
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
%
toc
K = [SIFT_Loc P3D(:,3)];
X = SIFT_Loc(:,1);
Z = P3D(:,3);
%%
LF1 = lighfield;
relx = [];
relz = [];
d=0;
GUI;
% figure(1)
% imshow(squeeze(LF1(6,6,:,:,1:3)));
% title('heres your image')
% n = input('enter the number of objects');

for i = 1:n
    figure(1)
    imshow(squeeze(LF1(6,6,:,:,1:3)));
    title(['select object number',i])
    rect = getrect;
    rect = [rect(1) rect(1)+rect(3)];
    temp = find(X<rect(2));
    temp1 = Z(temp);
    temp2 = X(temp);
    temp1 = temp1(find(temp2>rect(1)));
    temp1 = temp1(find(temp1>0));
    temp1 = temp1(find(temp1<0.19*2));
    z = median(temp1);
    relx = [relx rect(1):rect(2)];
    relz = [relz z*ones(1,ceil(rect(2)-rect(1))+1)];
end
close all
[~,~,x,y,~] = size(LF1);
[relx,I]=sort(relx,'ascend');
relz = relz(I);
mini = min(relx);
maxi = max(relx);
relx = [1:mini-1 relx maxi+1:x];
relz = [relz(1)*ones(1,ceil(mini-1)) relz relz(end)*ones(1,ceil(x-maxi))];
H1 = [0.0004 0 -0.072 ; 0 0.0004 -0.072 ; 0 0 1];
A = H1*[relx ;ones(size(relx));ones(size(relx))];
relu = A(1,:);
%%

H11 = inv(H1);
H2 = [0.0002 0 -0.0013 ; 0 0.0002 -0.0013 ; 0 0 1];
%%%%%
[l,k] = ndgrid(1:x,1:y);
k = reshape(k,[1,x*y]);
l = reshape(l,[1,x*y]);
%%%%%
[i,j] = ndgrid(1:11,1:11);
i = reshape(i,[1,11^2]);
j = reshape(j,[1,11^2]);
%%%%%
A = H1*[k;l;ones(1,x*y)];
u = A(1,:);
v = A(2,:);
%%%%%
A = H2*[i;j;ones(1,121)];
s = A(1,:);
t = A(2,:);
%%%%%
up = u;
vp = v;
%%
temp = reshape(u,[x y]);
temp = temp(1,:);
ztemp = interpn(relu,relz,temp,'spline',0);
%%
zp = reshape(ones(x,1)*ztemp,[1,x*y]);
%%%%%
F=0.19;
LF1 = im2double(LF1);
refocused=zeros(x,y,3);
    for i = 3:9
        for j = 3:9
            for k=1:3
                CurSlice=squeeze(LF1(j,i,:,:,k));
                u2 = up+s(i)*(1-F./zp);
                v2 = vp+s(j)*(1-F./zp);
                A = H11*[u2;v2;ones(1,x*y)];
                k2 = reshape(A(1,:),[x,y]);
                l2 = reshape(A(2,:),[x,y]);
                CurSlice = interpn(CurSlice,l2,k2,'spline',0);
                refocused(:,:,k)=CurSlice+refocused(:,:,k);
            end
        end
    end   
refocused = refocused/49;

figure(1)
imshow(refocused)
title('final image')

file = load('pen.mat');
image = file.LF;
LF1 = image(:,:,1:375,1:375,1:3);
% j = 6;
% epipolar = [];
% 
% for i = 1:11
%     temp = squeeze(image(6,i,:,:,1:3));
%     epipolar = cat(1,epipolar,temp(170,:,:));
% end

%%

[~,~,x,y,~] = size(LF1);
H1 = [0.0004 0 -0.072 ; 0 0.0004 -0.072 ; 0 0 1];
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
q1 = [-0.0388,0.0024,0.0753];
q2 = [0.0260,0.012,0.19];
T=[(q1(1)+q2(1))/2;0;(q1(3)+q2(3))/2];
zp = atan((400/0.78)*u)*(abs(q1(3)-q2(3)))/3.14;
T = T*ones(1,375^2);
A=[u;v;zp]+T;
up = A(1,:);
vp = A(2,:);
zp = A(3,:);
%%%%%
F=0.19;
LF1 = im2double(LF1);
refocused=zeros(x,y,3);
temp = zeros(x,y,3);
epipolar = [];
    for i = 1:11
        for j = 6
            for k=1:3
                CurSlice=squeeze(LF1(j,i,:,:,k));
                u2 = up+s(i)*(1-F./zp);
                v2 = vp+s(j)*(1-F./zp);
                A = H11*[u2;v2;ones(1,x*y)];
                k2 = reshape(A(1,:),[x,y]);
                l2 = reshape(A(2,:),[x,y]);
                CurSlice = interpn(CurSlice,l2,k2,'spline',0);
                temp(:,:,k) = CurSlice;
                refocused(:,:,k)=CurSlice+refocused(:,:,k);
            end
            epipolar = cat(1,epipolar,temp(170,:,:));
        end
    end   
refocused = refocused/49;
%%
figure(1)
imshow(epipolar)

%%
[~,~,x,y,~] = size(LF1);
H1 = [0.0004 0 -0.072 ; 0 0.0004 -0.072 ; 0 0 1];
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
Thetax = 0;
Thetay = -60;
Thetaz = 0;
T = [0.1561;0;0.035];
F = 0.19;
R = RotationMat(Thetax,Thetay,Thetaz);
T = T*ones(1,x*y);
A = R*[u;v;F*ones(1,x*y)]+T;
up = A(1,:);
vp = A(2,:);
zp = A(3,:);
%%%%%
LF1 = im2double(LF1);
refocused=zeros(x,y,3);
temp = zeros(x,y,3);
epipolar = [];

    for i = 1:11
        for j = 6
            for k=1:3
                CurSlice=squeeze(LF1(j,i,:,:,k));
                u2 = up+s(i)*(1-F./zp);
                v2 = vp+s(j)*(1-F./zp);
%                 u2 = (F./zp).*up+s(i)*(1-F./zp);
%                 v2 = vp+t(j)*(zp/F-1);
                A = H11*[u2;v2;ones(1,x*y)];
                k2 = reshape(A(1,:),[x,y]);
                l2 = reshape(A(2,:),[x,y]);
                CurSlice = interpn(CurSlice,l2,k2,'spline',0);
                temp(:,:,k) = CurSlice;
                refocused(:,:,k)=CurSlice+refocused(:,:,k);
            end
            epipolar = cat(1,epipolar,temp(170,:,:));
        end
    end   
refocused = refocused/49;

figure(2)
imshow(epipolar)

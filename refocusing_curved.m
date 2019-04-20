function [ refocused ] = refocusing_curved( LF )
[~,~,x,y,~] = size(LF);
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
LF = im2double(LF);
refocused=zeros(x,y,3);
temp = zeros(x,y,3);
epipolar = [];
    for i = 3:9
        for j = 3:9
            for k=1:3
                CurSlice=squeeze(LF(j,i,:,:,k));
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
end


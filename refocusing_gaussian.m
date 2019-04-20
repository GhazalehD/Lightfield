function [ refocused ] = refocusing_gaussian( LF )
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

mu=0.0016;
sigma=0.034/5;
pd = makedist('Normal',mu,sigma);
p = pdf(pd,u);
zp = p*0.21/max(p);
T = [0;0;0.14];
T = T*ones(1,375^2);
A=[u;v;zp]+T;
up = A(1,:);
vp = A(2,:);
zp = A(3,:);
%%%%%
F=0.19;
LF = im2double(LF);
refocused=zeros(x,y,3);
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
                refocused(:,:,k)=CurSlice+refocused(:,:,k);
            end
        end
    end   
refocused = refocused/49;

end


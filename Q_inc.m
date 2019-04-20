function [ A ] = Q_inc( image )

image = image(:,:,1:375,1:375,:);

[l,k] = ndgrid(1:2*375,1:2*375);
l2 = l/2;
k2 = k/2;

image = im2double(image);
A = zeros(11,11,375*2,375*2,3);
for i = 3:9
        for j = 3:9
            for k=1:3
                CurSlice=squeeze(image(j,i,:,:,k));
                A(j,i,:,:,k) = interpn(CurSlice,l2,k2,'spline',0);
            end
        end
end

end


function Refouced = refocusing (LF,pixel_shift)


[a,b,c,d,~]=size(LF);
centerx=(b+1)/2;
centery=(a+1)/2;
[yy,xx] = ndgrid(1:c,1:d);
LF = im2double(LF);
Refouced=zeros(c,d,3);
% h = 1;
% for pixel_shift = -1:0.2:1
    for i = 3:9
        for j = 3:9
            for k=1:3
                CurSlice=squeeze(LF(i,j,:,:,k));
                CurSlice = interpn(CurSlice, yy+(i-centery)*pixel_shift, xx+(j-centerx)*pixel_shift, 'linear', 0);
                Refouced(:,:,k)=CurSlice+Refouced(:,:,k);
            end
        end
    end
%     F(h)=im2frame(Refouced/max(max(max(Refouced))));
%     h=h+1;
% end

Refouced = Refouced/max(max(max(Refouced)));

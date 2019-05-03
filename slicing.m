
function [rgbimage]=slicing(Data3,a,f)
dim = size(Data3)./15;
slice1 = zeros(dim(1),dim(2),3);
% A = slice1;
% x = (1-a).*(1:0.1:433).*15+a.*(1:0.1:433);
% y = (1-a).*(1:0.1:625).*15+a.*(1:0.1:625);
% fakex=[1:433];
% fakey=[1:625];
for o = 1:3
%     A = interp2(fakex,fakey,Data3(:,:,o),x,y);
for kx = 1:dim(1)
    for ky = 1:dim(2)
        x = (1-a)*kx*15+a*kx;
        y = (1-a)*ky*15+a*ky;
%         I1 = Data3((ceil(x)), (ceil(y)), o);
%         w1 = ((ceil(x)-x)+(ceil(y)-y))/2;
        if (ceil(x)) == 0
            x=ceil(x)+1;
        end
        if (ceil(y)) == 0
            y=ceil(y)+1;
        end
%         I0 = Data3((ceil(x)-1), (ceil(y)-1), o);
%         w0 = ((x-(ceil(x)-1))+(y-(ceil(y)-1)))/2;
%         slice1(kx,ky,o) = 1/(f^2)*((I1*w1)+(I0*w0));
        slice1(kx,ky,o) = 1/(f^2)*Data3((ceil(x)), (ceil(y)), o);
        %slice1(kx,ky,o) = 1/(f^2)*A(x*10,y*10,o);        
    end
end
end

result1=ifft2(slice1(:,:,1));
result2=ifft2(slice1(:,:,2));
result3=ifft2(slice1(:,:,3));

result1 = (result1-min(result1(:)))/(max(result1(:))-min(result1(:)));
result2 = (result2-min(result2(:)))/(max(result2(:))-min(result2(:)));
result3 = (result3-min(result3(:)))/(max(result3(:))-min(result3(:)));

% R1 = interp2((result1-min(result1(:)))/(max(result1(:))-min(result1(:))),2);
% R2 = interp2((result2-min(result2(:)))/(max(result2(:))-min(result2(:))),2);
% R3 = interp2((result3-min(result3(:)))/(max(result3(:))-min(result3(:))),2);

rgbimage=mat2gray(abs(cat(3,result1,result2,result3)));

% rgbimage=abs(cat(3,result1,result2,result3));
%figure;imshow(abs(rgbimage))



end
clc; clear all; close all
%%
tic
%cd 'C:\Users\Alex\Desktop\Carney Research\LF TOOLBOX matlab\LFToolbox0.3_Samples1\LFToolbox0.3_Samples1\Images\Illum';
A=importdata('LorikeetHiding__Decoded.mat');
%% Process LF
LF = A.LF;
% NC = size(LF,5);
% NW = LF(:,:,:,:,end);
% for i = 1:NC
%     LF(:,:,:,:,i) = LF(:,:,:,:,i)./NW;
% end
% LF = max(0,LF);
% LF = min(1,LF);

LF2 = A.LF(:,:,:,:,1);% using one of the color

LF_size = size(LF2);
%%
% a = 0; % scaling factor
%B = [a 0 1-a 0; 0 a 0 1-a; 0 0 1 0; 0 0 0 1]; % shearing
%% 4D fourier transform
% x = 1:LF_size(1);
% y = 1:LF_size(2);
% u = 1:LF_size(3);
% v = 1:LF_size(4);
% LF3 =zeros(size(LF2));

%% 4D transform
% for i = x
%     AA(i,:,:,:,1:3)  = fft(LF(i,:,:,:,1:3));
% end
% 
% BB=zeros(size(AA));
% CC=BB;
% DD=CC;
% E=DD;

% with fft shift
% for j = y
%     BB(:,j,:,:,1:3)  = fftshift(fft(AA(:,j,:,:,1:3)));
% end
% 
% for k = u
%     CC(:,:,k,:,1:3)  = fftshift(fft(BB(:,:,k,:,1:3)));
% end
% 
% for l = v
%     DD(:,:,:,l,1:3)  = fftshift(fft(CC(:,:,:,l,1:3)));
% end

% without fftshift
% for j = y
%     BB(:,j,:,:,1:3)  = fft(AA(:,j,:,:,1:3));
% end
% 
% for k = u
%     CC(:,:,k,:,1:3)  = fft(BB(:,:,k,:,1:3));
% end
% 
% for l = v
%     DD(:,:,:,l,1:3)  = fft(CC(:,:,:,l,1:3));
% end
% 
% E = fftshift(DD);
%%
 %figure;imshow(squeeze(E2(:,:,:,:,1:3)));
 %figure;imshow(squeeze(E(:,:,250,250,1:3)));
%%
% tic
% B2 = a^2.*[1 0 1-1/a 0; 0 1 0 1-1/a; 0 0 1/a 0; 0 0 0 1/a];
% for l = v 
%     for k = u
%         for j = y
%             for i = x
%                 z = B2*[i;j;k;l];
%                 aaa = z(1);
%                 bbb = z(2);
%                 ccc = z(3);
%                 ddd = z(4);
%                 for ggg = 1:3
%                     F(aaa,bbb,ccc,ddd,ggg) = E(i,j,k,l,ggg);
%                 end
%             end
%         end
%     end
% end
% toc

%% Light field
% lenslet image
dim1 =433;
dim2 =625;
Data2 = zeros([15*dim1 15*dim2 3]); 
Data3 = Data2; 
for i = 1:dim1
    for j = 1:dim2
        Data{i,j} = squeeze(LF(:,:,i,j,1:3));
        Data4((i+15*(i-1)-(i-1)):15*i,((j+15*(j-1)-(j-1)):15*j),:) = Data{i,j}; 
        Data2((i+15*(i-1)-(i-1)):15*i,((j+15*(j-1)-(j-1)):15*j),:) = fftshift(fft2(Data{i,j}));
    end
end

%%
dim1 =434;
dim2 =625;
% subaperture image

for i = 1:15
    for j = 1:15
        Data5{i,j} = squeeze(LF(i,j,:,:,1:3));
        Data6((i+dim1*(i-1)-(i-1)):dim1*i,((j+dim2*(j-1)-(j-1)):dim2*j),:)=Data5{i,j}; 
        Data7((i+15*(i-1)-(i-1)):15*i,((j+15*(j-1)-(j-1)):15*j),:) = fftshift(fft2(Data{i,j}));
    end
end
%% plot three different graphs

% first color
Data3(:,:,1) = (fft2(squeeze(Data4(:,:,1))));% after 2D fft on each 15*15 then stacked.

% figure;imshow(squeeze(Data3(:,:,1)),[0 2^30]) % fourier transform in total
% figure;imshow(squeeze(Data2(:,:,1:3)),[0 2^16]) % fourier transform in lenslet image
% figure;imshow(squeeze(Data4(:,:,1:3)),[0 2^16]) % original lenslet image

% second color

Data3(:,:,2) = (fft2(squeeze(Data4(:,:,2))));% after 2D fft on each 15*15 then stacked.

% figure;imshow(squeeze(Data3(:,:,2)),[0 2^30]) % fourier transform in total

% third color 

Data3(:,:,3) = (fft2(squeeze(Data4(:,:,3))));% after 2D fft on each 15*15 then stacked.

% figure;imshow(squeeze(Data3(:,:,1)),[0 2^30]) % fourier transform in total

%% Operating FFT on subaperture image
Data33(:,:,1) = (fft2(squeeze(Data6(:,:,1))));% after 2D fft on each 15*15 then stacked.
Data33(:,:,2) = (fft2(squeeze(Data6(:,:,2))));% after 2D fft on each 15*15 then stacked.
Data33(:,:,3) = (fft2(squeeze(Data6(:,:,3))));% after 2D fft on each 15*15 then stacked.

%% fourier slice

% Slice1=Data3((6510/2-433/2):(6510/2+433/2),(9374/2-624/2):(9374/2+624/2),1:3);
% %figure;imshow(squeeze(Slice1(:,:,1:3)))
% result1=ifft2(squeeze(Slice1(:,:,1:3)));
% figure;imshow(squeeze(result1(:,:,1:3)))

% a = 1.3;
% f = 1;
% b = 14;
% 
% slice1 = zeros(433,625,3);
% for o = 1:3
% for kx = 1:433
%     for ky = 1:625
%         slice1(kx,ky,o) = 1/(f^2)*Data3(round(abs((1-a)*kx*15+a*kx)), round(abs((1-a)*ky*15+a*ky)), o);
%     end
% end
% end 
% figure(2)
% result1=ifft2(slice1(:,:,1));
% imshow(abs(result1), [0 2^b])
% 
% figure(3)
% result2=ifft2(slice1(:,:,2));
% imshow(abs(result2), [0 2^b])
% 
% figure(4)
% result3=ifft2(slice1(:,:,3));
% imshow(abs(result3), [0 2^b])
% 
% % figure
% % result1=ifft2(squeeze(slice1(:,:,1:3)));
% % imshow(abs(squeeze(result1(:,:,1:3))))
% R1 = (result1-min(result1(:)))/(max(result1(:))-min(result1(:)));
% R2 = (result2-min(result2(:)))/(max(result2(:))-min(result2(:)));
% R3 = (result3-min(result3(:)))/(max(result3(:))-min(result3(:)));
% 
% rgbimage=cat(3,R1,R2,R3);
% % for i = 1:3
% %     R4(:,i) = (rgbimage(:,i)-min(rgbimage(:,i)))/(max(rgbimage(:,i))-min(rgbimage(:,i)));
% % end
% 
% figure;imshow(abs(rgbimage))


%% create a video of the refocus images 
num = 10; % number of images
alpha = linspace(0.9,1,num); % alpha for refocus
f=1; % alpha = f'/f, which is the fraction 
image=cell(num,1); % the image series
position = [50 50]; % position of the text box 


% create text marker for each image
text_str = cell(num,1);
for ii=1:num
   text_str{ii} = ['alpha = ' num2str(alpha(ii),'%0.4f')];
end

% insert text to each image
for i = 1:num
    image{i} = slicing(Data3,alpha(i),f);
    %image{i} = slicing(Data33,alpha(i),f);
    image{i} = insertText(image{i}, position, text_str{i},'FontSize',18,'TextColor','white');
end


 % create the video writer with 1 fps
 writerObj = VideoWriter('report1.avi');
 writerObj.FrameRate = 1;
 % set the seconds per image
%  secsPerImage = [5 10 15];
 % open the video writer
 open(writerObj);
 % write the frames to the video
 for u=1:length(image)
     % convert the image to a frame
     % frame = im2frame(image{u});
     
         writeVideo(writerObj, image{u});
 
 end
 % close the writer object
 close(writerObj);
 
 toc
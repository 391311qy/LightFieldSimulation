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

%% Apply eq 4.2
alpha = 1;
tic
AA = zeros(15,15);
for k = 1:434
    for l = 1:625
                x = [1:15].*(1-1/alpha)+k/alpha;
                y = [1:15].*(1-1/alpha)+l/alpha;
                AA(:,:) = Data6((ceil(x)),ceil((y)),1);
        
    end
end

%%
dim1 =434;
dim2 =625;

 [XX YY] = meshgrid(1:dim2, 1:dim1);
 
%% epipolar images
% for i = 1:434
%     for j =1:15
%     AA{i,j} = squeeze(LF(:,j,:,i,1:3));
%     end
% end
% 
% for j = 1:15
% BB{j} = vertcat(AA{:,j});
% end
% 
% BB = horzcat(BB{:});
% imshow(BB)
% %%
% SS = uint16(zeros(15,434,3));
% for i =1:434
% SS = SS + (AA{i,5})./434;
% end
% 
%%

Data3(:,:,1) = (fft2(squeeze(Data4(:,:,1))));% after 2D fft on each 15*15 then stacked.
Data3(:,:,2) = (fft2(squeeze(Data4(:,:,2))));% after 2D fft on each 15*15 then stacked.
Data3(:,:,3) = (fft2(squeeze(Data4(:,:,3))));% after 2D fft on each 15*15 then stacked.

%% other people algorithm
num = 20;
alpha = linspace(-1,1,num); % alpha for refocus
position = [50 50]; % position of the text box 

% create text marker for each image
text_str = cell(num,1);
for ii=1:num
   text_str{ii} = ['alpha = ' num2str(alpha(ii),'%0.4f')];
end

% insert text to each image
for i = 1:num
    image{i} = refocusLightField(double(A.LF),alpha(i));
    image{i} = insertText(image{i}, position, text_str{i},'FontSize',18,'TextColor','white');
end

 % create the video writer with 1 fps
 writerObj = VideoWriter('copycode.avi');
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

 
 
 %% write image original code
% num = 10; % number of images
% alpha = linspace(0.97,1,num); % alpha for refocus
% f=1; % alpha = f'/f, which is the fraction 
% image=cell(num,1); % the image series
% position = [50 50]; % position of the text box 
% 
% 
% % create text marker for each image
% text_str = cell(num,1);
% for ii=1:num
%    text_str{ii} = ['alpha = ' num2str(alpha(ii),'%0.4f')];
% end
% 
% % insert text to each image
% for i = 1:num
%     image{i} = slicing(Data3,alpha(i),f);
%     %image{i} = slicing(Data33,alpha(i),f);
%     image{i} = insertText(image{i}, position, text_str{i},'FontSize',18,'TextColor','white');
% end
% 
% 
%  % create the video writer with 1 fps
%  writerObj = VideoWriter('report1.avi');
%  writerObj.FrameRate = 1;
%  % set the seconds per image
% %  secsPerImage = [5 10 15];
%  % open the video writer
%  open(writerObj);
%  % write the frames to the video
%  for u=1:length(image)
%      % convert the image to a frame
%      % frame = im2frame(image{u});
%      
%          writeVideo(writerObj, image{u});
%  
%  end
%  % close the writer object
%  close(writerObj);
%  
 toc
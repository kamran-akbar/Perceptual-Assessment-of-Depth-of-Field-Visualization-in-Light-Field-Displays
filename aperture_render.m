clear all;
clc;

tic

dirInStr = string(zeros(1, 3));
dirOutStr = string(zeros(1, 3));
for k = 1:6
    dirInStr(k) = strcat("D:/Scenes/Scene ", num2str(k), "/FP/"); 
    dirOutStr(k) = strcat("D:/Scenes/Scene ", num2str(k), "/AF/");
end


cam.colRes = 1280;
cam.rowRes = 720;
cam.nrCamHor = 699;  % [0000:0660]
cam.nrCamVer = 40;  % [0000:0175]

aper.ranges = 0:3:15;   % radius ... 33 --> 13 GBy memory
verCenter = floor(cam.nrCamVer/2);


%% pre-calculate image indexes for each apperture (in terms of difference
% with the smaller apperture)
dir = 5;
%for dir = 3:4 %squeeze(size(dirInStr, 2))
    [xGrid, yGrid] = meshgrid(floor(-max(aper.ranges)) : ceil(max(aper.ranges)));
    xyRadius = sqrt(xGrid(:).^2 + yGrid(:).^2);
    
    aper.xInd{length(aper.ranges)} = [];
    aper.yInd{length(aper.ranges)} = [];
    xAll = []; yAll = [];
    
    aper.xInd{1} = xGrid(xyRadius <= aper.ranges(1));
    aper.yInd{1} = yGrid(xyRadius <= aper.ranges(1));
    aper.nrPoints(1) = length(aper.xInd{1});
    
    for k = 2:length(aper.ranges)
          rangeAmount = aper.ranges(k);
          
          xIn = xGrid((xyRadius <= rangeAmount) & (xyRadius > aper.ranges(k-1)));
          yIn = yGrid((xyRadius <= rangeAmount) & (xyRadius > aper.ranges(k-1)));
    
          aper.xInd{k} = xIn;
          aper.yInd{k} = yIn;
          
          aper.nrPoints(k) = aper.nrPoints(k-1) + length(xIn);
    end
    %% preload images for position 1 --> assumes that apperture is applied to 
    % all cameras are caluclated starting with camera 0 
    imIn = zeros(2*max(aper.ranges)+1, 2*max(aper.ranges)+1, cam.rowRes, cam.colRes, 3, 'uint8');
    imInInd = NaN * zeros(2*max(aper.ranges)+1, 1);
    
    for m = 1:max(aper.ranges)+1                    % x
       for k = -max(aper.ranges):max(aper.ranges)   % y
          imStr = strcat(dirInStr(dir), 'img_V', num2str(verCenter + k, '%04d'), '_H', num2str(m-1, '%04d'), '.png');
           if exist(imStr, 'file')
              imTmp = imread(imStr);
              if (isa(imTmp, "uint16"))
                 imTmp = uint8(double(imTmp) * (2^8-1)/(2^16-1));
              end
              imIn(m+max(aper.ranges), k+max(aper.ranges)+1, :, :, :) = imTmp;          
           end
       end
    end
    imInInd(max(aper.ranges)+1 : end) = 0;
    %% average images over appertures
    xyCenter = max(aper.ranges)+1;
    % all camera positions
    for m = 1:cam.nrCamHor
       disp([m cam.nrCamHor]);
     
        im = zeros(cam.rowRes, cam.colRes, 3, 'double');
        counterCorrection = 0;   % compensating for missing images on the edges
       
        % all appertures
        for k = 1:length(aper.ranges)  
          
           % all images
           for n = 1:length(aper.xInd{k})
              if isnan(imInInd(xyCenter+aper.xInd{k}(n)))
                 counterCorrection = counterCorrection + 1;
              else
                imTmp = squeeze(imIn(xyCenter + aper.xInd{k}(n), xyCenter + aper.yInd{k}(n), :, :, :));
                im = im + double(imTmp);
              end
           end
          
          imOut = uint8(im / (aper.nrPoints(k) - counterCorrection));
%           imOutText = insertText(imOut, [500 50], num2str(aper.ranges(k),'%02d'), 'TextColor', 'w', 'FontSize', 30, 'BoxColor', [0 0 0]);
            if dir == 1
                imwrite(imOut, strcat(dirOutStr(dir), 'A', num2str(k-1, '%04d'), '_I', num2str(cam.nrCamHor-m, '%04d'), '.png'));
            else
                imwrite(imOut, strcat(dirOutStr(dir), 'A', num2str(k-1, '%04d'), '_I', num2str(m-1, '%04d'), '.png'));
            end
        end
      
        % prepare data for next camera position (add one more column, shift
        % everything to left)
        imIn = circshift(imIn, -1);
        imInInd = circshift(imInInd, -1);
        imInInd(end) = NaN;
        
        for k = -max(aper.ranges):max(aper.ranges)   % y
          imStr = strcat(dirInStr(dir), 'img_V', num2str(verCenter + k, '%04d'), '_H', num2str(m + max(aper.ranges), '%04d'), '.png'); % m is effectively next image (m+1)
           if exist(imStr, 'file')
              imTmp = imread(imStr);
              if (isa(imTmp, "uint16"))
                 imTmp = uint8(double(imTmp) * (2^8-1)/(2^16-1));
              end
              imIn(end, k+max(aper.ranges)+1, :, :, :) = imTmp;
              imInInd(end) = 0;
           end
        end
       
        toc
        
     end
%end
toc
% data = load_images(dirInStr, 'img_V', [0 1 cam.nrCamVer-1 4], '_H0000.png', 'uint8');

% write_images('Vessels_orig331', dataCam(1:2:end,:,:,:));

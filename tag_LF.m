%%
for k = 1:5
    for im_num = 0:698
        I = imread(strcat("D:/Scenes/Scene 5/AF_old/", 'A', num2str(k, '%04d'), '_I', num2str(im_num, '%04d'), '.png'));
        imwrite(I, strcat("D:/Scenes/Scene 5/AF/", 'A', num2str(k*3, '%04d'), '_I', num2str(im_num, '%04d'), '.png'));
    end
end

%%
clear all;
clc;


dirInStr = string(zeros(1, 5));
dirOutStr = string(zeros(1, 5));

for k = 5
    dirInStr(k) = strcat("D:/Scenes/Scene ", num2str(k), "/AF/"); 
    dirOutStr(k) = strcat("D:/Scenes/Scene ", num2str(k), "/AF_tagged/");
end
for aper = 0:3:15
    for k = 5
        filename =strcat('A', num2str(aper, '%04d'), '_I', num2str(19, '%04d'), '.png');
        I = imread(strcat(dirInStr(k), filename));
        imOutText = insertText(I, [20 650], num2str(aper,'%02d'), 'TextColor', 'w', 'FontSize', 15, 'BoxColor', [0 0 0]);
        imshow(imOutText)
        imwrite(imOutText, strcat(dirOutStr(k), filename));
    end
end
%%
clear all;
clc;


dirInStr = string(zeros(1, 3));
dirOutStr = string(zeros(1, 3));

for k = 1:3
    dirInStr(k) = strcat("D:/Scenes/Scene ", num2str(k), "/HPO/"); 
    dirOutStr(k) = strcat("D:/Scenes/Scene ", num2str(k), "/HPO_tagged/");
end
for aper = 0:19
    for k = 1
        filename = strcat('order_', num2str(aper, '%04d'), 'I_', num2str(680, '%04d'), '.png');
        I = imread(strcat(dirInStr(k), filename));
        imOutText = insertText(I, [20 650], num2str(aper,'%02d'), 'TextColor', 'w', 'FontSize', 15, 'BoxColor', [0 0 0]);
        imshow(imOutText)
        filename2 = strcat('order_', num2str(aper, '%04d'), '_I', num2str(20, '%04d'), '.png');
        imwrite(imOutText, strcat(dirOutStr(k), filename2));
    end
end
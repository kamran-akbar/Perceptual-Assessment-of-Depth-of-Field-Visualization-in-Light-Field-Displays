%%
clc;
clear all;
close all;
%%
aperture_scene = [9, 9, 0, 9, 9, 9];
% metric_scene = [0.647, 0.15, 3.18, 0.13, 0.42, 0.25]; %w1w2w3
% metric_scene = [0.136, 0.075, 0.388, 0.0935, 0.2642, 0.063]; %w1w2w3_0.2DoF
metric_scene = [0.19, 0.094, 0.52, 0.093, 0.22, 0.067]; %w1w2w3_0.1DoF
% metric_scene = [1.18, 0.4, 2.5, 0.24, 0.54, 0.41];    %w1(w2+w3)
% metric_scene = [1.8, 1.49, 4.35, 0.96, 1.95, 1.52];   %w1+w2+w3
labels = ["V" "R" "D" "T" "Z" "L"];
aperture_scene_val = [1 1 1];
metric_scene_val = [0.489, 0.02, 0.094];
rho = corr(aperture_scene',metric_scene');

modelFun = @(b,x) b(3)./(1+exp(b(1).*(x-b(2))));
start = [48 0.4 9];
xx = 0.0:0.000001:0.7;
% y = start(3)./(1+exp(start(1).*(xx-start(2))));
% figure;
% plot(xx, y);

% figure;
% nlm = fitnlm(metric_scene,aperture_scene1,modelFun,start);
% scatter(metric_scene, aperture_scene1);
% hold on;
% line(xx,predict(nlm,xx'),'linestyle','--','color','k');
% title("Dining blur, Lab sharp");
% xlabel("Metric");
% ylabel("Aperture");
% hold off;

% figure;
% nlm_blur = fitnlm(metric_scene,aperture_scene2,modelFun,start);
% scatter(metric_scene, aperture_scene2);
% hold on;
% line(xx,predict(nlm_blur,xx'),'linestyle','--','color','k');
% title("Dining blur, Lab blur");
% xlabel("Metric");
% ylabel("Aperture");
% hold off;

figure;
nlm = fitnlm(metric_scene,aperture_scene,modelFun,start);
scatter(metric_scene, aperture_scene);
hold on;
line(xx,predict(nlm,xx'),'linestyle','--','color','k');
title("Dining sharp, Lab sharp");
xlabel("Metric");
ylabel("Aperture");
hold off;
%%
min_pref_aperture =min(aperture_scene);
max_pref_aperture =max(aperture_scene);
figure;
blur_model =predict(nlm,xx');
line(xx, blur_model,'linestyle','--','color','b');
hold on;
sharp_model = predict(nlm,xx');
line(xx,sharp_model,'linestyle','--','color','b');
xlabel("$f$", 'Interpreter','latex', 'FontSize', 14);
ylabel("Aperture Radii");
hold on;
points = zeros(max_pref_aperture-min_pref_aperture, 2);
points(1, :) = [0 9]'; 
count = 2;
% for aper =min_pref_aperture+1:max_pref_aperture-1
%     a = xx(abs(sharp_model-aper) <0.1);
%     b = xx(abs(blur_model-aper) <0.1);
%     mid = (a(ceil(size(a,2) / 2)) + b(ceil(size(b,2) / 2))) / 2;
%     points(count, 1) = mid;
%     points(count, 2) = aper;
%     count = count + 1;
% end
% points(count, :) = [1 0]';
% scatter(points(:, 1), points(:, 2), 'color', 'g');
% hold on;
% start(3) = nlm_sharp.Coefficients.Estimate(3); 
hold on;
scatter(metric_scene, aperture_scene, 'color', 'r');
textPosX = metric_scene -0.005;
textPosX(2) = textPosX(2) - 0.01;
textPosX(4) = textPosX(4) + 0.01;
textPosY = aperture_scene + 0.7;
text(textPosX, textPosY, labels, "FontSize", 10);
ylim([0 15])
hold off;
%%
aperture_scene_val = [0 12 9];
metric_scene_val = [0.55, 0.12, 0.041];
labels_val = ["C" "G" "F"];
figure;
blur_model =predict(nlm,xx');
line(xx, blur_model,'linestyle','--','color','b');
hold on;
scatter(metric_scene_val, aperture_scene_val, 'color', 'b');
text(metric_scene_val, aperture_scene_val+0.5,labels_val)
xlabel("$f$", 'Interpreter','latex', 'FontSize', 14);
ylabel("Aperture Radii");
xlim([0 0.7]);
ylim([0 15])
hold off;
v_estimate = aperture_scene_val;
v_true = predict(nlm, metric_scene_val')';
eps = sum(abs(v_estimate - v_true))/length(v_true);
% eps = sum(abs(floor((v_estimate - v_true) ./ v_true)))/length(v_true);



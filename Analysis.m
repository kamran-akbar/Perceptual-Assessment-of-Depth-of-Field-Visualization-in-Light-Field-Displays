
%% Validation study
clc;
close all;
clear all;
participantNumber = 14;
SceneNumber = 3;
comparisonPair = 9;
outlier_indices = [];
validParticipant = participantNumber - size(outlier_indices, 2);
scene_flower = zeros(validParticipant, comparisonPair, 2);
scene_camper = zeros(validParticipant, comparisonPair, 2);
scene_garden = zeros(validParticipant, comparisonPair, 2);
counter = 1;
for ii = 1:participantNumber
    if (~ismember(ii, outlier_indices))
        load(strcat('./playlists_evaluation/randomKeyGen/', 'G', num2str(ii), 'S7.mat'));
        scene_flower(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists_evaluation/randomKeyGen/', 'G', num2str(ii), 'S8.mat'));
        scene_camper(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists_evaluation/randomKeyGen/', 'G', num2str(ii), 'S9.mat'));
        scene_garden(counter, :, :) = squeeze(randomized_pair);
        counter = counter + 1;
    end
end
[num,txt,raw] = xlsread('Validation Study.xlsx') ;
people = raw;
people = people(2:end, 2:end);
people = people(1:participantNumber, :);
people(outlier_indices, :) = [];
load('sceneRandomEvaluation.mat');
randScene = randSceneEve(1:participantNumber, :);
people_order = cell(size(people));
randScene(outlier_indices, :) = [];
for p = 1:size(randScene, 1)
    for s = 1:size(randScene, 2)
        s_idx = randScene(p, s);
        range_rnd = (comparisonPair * (s_idx - 7) + 1):(comparisonPair * (s_idx - 6));
        range = (comparisonPair * (s - 1) + 1):(comparisonPair * s);
        people_order(p, range_rnd) = people(p, range);
    end
end
max_aper = 15;
min_aper = 0;
interval = 3;
aperture = min_aper:interval:max_aper;

people_order_ = people_order(~cellfun('isempty', people_order));
people_order_ = reshape(people_order_, participantNumber, comparisonPair*size(randScene, 2));
W_flower = computeRankings(people_order_, scene_flower, 1:comparisonPair, interval);
W_scenes(:, :, 1) = W_flower;
W_camper = computeRankings(people_order_, scene_camper, comparisonPair+1:comparisonPair*2, interval);
W_scenes(:, :, 2) = W_camper;
W_garden = computeRankings(people_order_, scene_garden, comparisonPair*2+1:comparisonPair*3, interval);
W_scenes(:, :, 3) = W_garden;
choice_blur = zeros(size(people, 1), 3);
probs = zeros(SceneNumber, 3);
names = ["Flower" "Camper" "Garden"];
scn_p_aper = zeros(size(names, 2), size(people, 1), 2);
scenes_random = zeros(size(people, 1), comparisonPair, 2, SceneNumber);
scenes_random(:, :, :, 1) = scene_flower;
scenes_random(:, :, :, 2) = scene_camper;
scenes_random(:, :, :, 3) = scene_garden;

for jj = 1:SceneNumber
     for personID = 1:size(people, 1)
        W = computeRankingsSingle(people_order, scenes_random(:, :, :, jj) , ...
        comparisonPair*(jj-1)+1:comparisonPair*jj, interval, personID);
        P = ones(1, size(aperture, 2));
        P = Bradley_Terry_old(W, P);
        [P_sorted, idx] = sort(P);
        choice_blur(personID, 1) = (idx(end)-1)*interval;
        choice_blur(personID, 2) = (idx(end-1)-1)*interval;
        choice_blur(personID, 3) = (idx(end-2)-1)*interval;
        scn_p_aper(jj, personID, 1) = choice_blur(personID, 1);
        scn_p_aper(jj, personID, 2) = choice_blur(personID, 2);
        scn_p_aper(jj, personID, 3) = choice_blur(personID, 3);
     end
    for ii = min_aper:interval:max_aper
        probs(ii/interval+1, 1) = numel(find(choice_blur(:, 1) == ii))/size(people, 1);
        probs(ii/interval+1, 2) = numel(find(choice_blur(:, 2) == ii))/size(people, 1);
        probs(ii/interval+1, 3) = numel(find(choice_blur(:, 3) == ii))/size(people, 1);
        probs(ii/interval+1, 4) = sum(probs(ii/interval+1, :));
    end
    figure;
    bar(min_aper:interval:max_aper, probs);
    legend('First choice', 'Second choice');
    title(names(jj));
    xlabel("Aperture Radii")
    ylabel("Probability")
end
figure;
subplot(2,1,1);
b=bar3(1:SceneNumber, squeeze(scn_p_aper(:, :, 1)));
title("Aperture preference for each person and scene");
xlabel("Participant Number")
ylabel("Scene Number")
colorbar;
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
std_aperture = std(scn_p_aper(:, :, 1));
mean_aperture = mean(scn_p_aper(:, :, 1));
I = zeros(16, size(people, 1));
for kk=1:size(people, 1)
    if std_aperture(1, kk) == 0
        std_aperture(1, kk) = std_aperture(1, kk) + 1;
    end
    I(:, kk) = normpdf((0:15)',mean_aperture(1, kk),std_aperture(1, kk));
end
[sorted_mean, idx] = sort(mean_aperture);
% [sorted_std, idx] = sort(std_aperture(idx));
I_sorted = I(:, idx);
subplot(2,1,2);
[X,Y] = meshgrid(1:size(people, 1), 0:max_aper);
plot3(X, Y, I_sorted);
colorbar;
title("Normal Distribution of Aperture Preference Among Scenes");
xlabel("Participant Number")
ylabel("Aperture")
% c = choice_blur(find(choice_blur(:, 1) == 0), :);
% c_prob = numel(find(c(:, 2) == 9))/numel(c(:, 1));
% p_0 = numel(find(choice_blur(:, 1) == 0 | choice_blur(:, 2) == 0)) / participantNumber;
% p_9 = numel(find(choice_blur(:, 1) == 9 | choice_blur(:, 2) == 9)) / participantNumber;
%%
selected_prob = zeros(SceneNumber, max_aper/interval+1);
log_Lambda = zeros(SceneNumber, max_aper/interval+1, max_aper/interval+1);
for jj=1:SceneNumber
    P = ones(1, size(aperture, 2));
    P = Bradley_Terry_Newman(W_scenes(:, :, jj), P);
    figure;
    stem(aperture, P);
    hold on;
    scatter(aperture, P, "blue", "filled")
    hold on;
    stem(aperture(1), P(1), "red", "filled");
    hold on;
    if jj ~= 2
        stem(aperture(4), P(4), "red", "filled");
    end
    hold off;
    xticks(aperture)
    title(names(jj)); 
    xlabel("Aperture Radii")
    ylabel("$q_i$", 'Interpreter', 'latex', 'FontSize', 17)
    ylim([0 9]);
    for r1=1:max_aper/interval+1
        for r2=1:max_aper/interval+1
            if r1~=r2
                prob1 = (P(r1))/((P(r1)) + (P(r2)));
                prob2 = 1-prob1;
                n12 = W_scenes(r1, r2, jj);
                n21 = W_scenes(r2, r1, jj);
                Lambda= 0.5^(n12+n21)/(prob1^n12 * prob2^n21);
                temp = log(0.5^(n12+n21)) - log(prob1^n12 * prob2^n21);
                 if(temp>0)
                    prob1 = n12/(n12+n21);
                    prob2 = 1-prob1;
                    temp = log(0.5^(n12+n21)) - log(prob1^n12 * prob2^n21);
                end
                log_Lambda(jj, r1, r2) = abs(2*temp);
            end
        end
    end
    
    tempLambda = squeeze(log_Lambda(jj, :, :));
    a = sum(log(P));
    for rr=1:6
        selected_prob(jj, rr) = sum(squeeze(W_scenes(rr, :, jj)))...
            /(size(people,1)*comparisonPair);
    end
end
%% Correlational study
close all;
clear all;
clc;
% loading random key for each scene
participantNumber = 31;
comparisonPair = 9;
SceneNumber =6;
outlier_indices = [5 7 12 28];
validParticipant = participantNumber - size(outlier_indices, 2);
scene_vessel = zeros(validParticipant, comparisonPair, 2);
scene_dragon = zeros(validParticipant, comparisonPair, 2);
scene_zoo = zeros(validParticipant, comparisonPair, 2);
scene_toy = zeros(validParticipant, comparisonPair, 2);
scene_lab = zeros(validParticipant, comparisonPair, 2);
scene_dining = zeros(validParticipant, comparisonPair, 2);
counter = 1;
for ii = 1:participantNumber
    if (~ismember(ii, outlier_indices))
        load(strcat('./playlists/randomGenKey/', 'O', num2str(ii), 'S1.mat'));
        scene_vessel(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists/randomGenKey/', 'O', num2str(ii), 'S2.mat'));
        scene_dragon(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists/randomGenKey/', 'O', num2str(ii), 'S3.mat'));
        scene_zoo(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists/randomGenKey/', 'O', num2str(ii), 'S4.mat'));
        scene_toy(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists/randomGenKey/', 'O', num2str(ii), 'S5.mat'));
        scene_lab(counter, :, :) = squeeze(randomized_pair);
        load(strcat('./playlists/randomGenKey/', 'O', num2str(ii), 'S6.mat'));
        scene_dining(counter, :, :) = squeeze(randomized_pair);
        counter = counter + 1;
    end
end

% read votes and ordered them
[num,txt,raw] = xlsread('Correlation Study.xlsx') ;
people = raw;
people = people(2:end, 2:end);
people = people(1:participantNumber, :);
people(outlier_indices, :) = [];
load('sceneRandom2.mat');
randScene = randScene(1:participantNumber, :);
poeple_order = cell(size(people));
randScene(outlier_indices, :) = [];
for p = 1:size(randScene, 1)
    for s = 1:size(randScene, 2)
        s_idx = randScene(p, s);
        range_rnd = (comparisonPair * (s_idx - 1) + 1):(comparisonPair * s_idx);
        range = (comparisonPair * (s - 1) + 1):(comparisonPair * s);
        poeple_order(p, range_rnd) = people(p, range);
    end
end
% compute ranks
max_aper = 15;
min_aper = 0;
interval = 3;
aperture = min_aper:interval:max_aper;
W_scenes = zeros(size(aperture, 2), size(aperture, 2), SceneNumber);
W_vessel = computeRankings(poeple_order, scene_vessel, 1:comparisonPair, interval);
W_scenes(:, :, 1) = W_vessel;
W_dragon = computeRankings(poeple_order, scene_dragon, comparisonPair+1:comparisonPair*2, interval);
W_scenes(:, :, 2) = W_dragon;
W_zoo = computeRankings(poeple_order, scene_zoo, comparisonPair*2+1:comparisonPair*3, interval);
W_scenes(:, :, 3) = W_zoo;
W_toy = computeRankings(poeple_order, scene_toy, comparisonPair*3+1:comparisonPair*4, interval);
W_scenes(:, :, 4) = W_toy;
W_lab = computeRankings(poeple_order, scene_lab, comparisonPair*4+1:comparisonPair*5, interval);
W_scenes(:, :, 5) = W_lab;
W_dining = computeRankings(poeple_order, scene_dining, comparisonPair*5+1:comparisonPair*6, interval);
W_scenes(:, :, 6) = W_dining;
choice_blur = zeros(size(people, 1), 3);
probs = zeros(SceneNumber, 2);
names = ["Vessel" "Dragon" "Zoo" "Toy" "Laboratory" "Dining Room"];
scn_p_aper = zeros(size(names, 2), size(people, 1), 2);
scenes_random = zeros(size(people, 1), comparisonPair, 2, SceneNumber);
scenes_random(:, :, :, 1) = scene_vessel;
scenes_random(:, :, :, 2) = scene_dragon;
scenes_random(:, :, :, 3) = scene_zoo;
scenes_random(:, :, :, 4) = scene_toy;
scenes_random(:, :, :, 5) = scene_lab;
scenes_random(:, :, :, 6) = scene_dining;
for jj = 1:SceneNumber
     P_total = zeros(size(people, 1), size(aperture, 2));
     aperture_total = zeros(size(people, 1), size(aperture, 2));
     for personID = 1:size(people, 1)
        W = computeRankingsSingle(poeple_order, scenes_random(:, :, :, jj) , ...
            comparisonPair*(jj-1)+1:comparisonPair*jj, interval, personID);
        P = ones(1, size(aperture, 2));
        P = Bradley_Terry_old(W, P);
        [P_sorted, idx] = sort(P);
        P_total(personID, :) = P_sorted;
        aperture_total(personID, :) = (idx-1)*interval;
        choice_blur(personID, 1) = (idx(end)-1)*interval;
        choice_blur(personID, 2) = (idx(end-1)-1)*interval;
        choice_blur(personID, 3) = (idx(end-2)-1)*interval;
        choice_blur(personID, 4) = (idx(end-3)-1)*interval;
        choice_blur(personID, 5) = (idx(end-4)-1)*interval;
        choice_blur(personID, 6) = (idx(end-5)-1)*interval;
        scn_p_aper(jj, personID, 1) = choice_blur(personID, 1);
        scn_p_aper(jj, personID, 2) = choice_blur(personID, 2);
        scn_p_aper(jj, personID, 3) = choice_blur(personID, 3);
     end
    for ii = min_aper:interval:max_aper
        probs(ii/interval+1, 1) = numel(find(choice_blur(:, 1) == ii))/size(people, 1);
        probs(ii/interval+1, 2) = numel(find(choice_blur(:, 2) == ii))/size(people, 1);
        probs(ii/interval+1, 3) = numel(find(choice_blur(:, 3) == ii))/size(people, 1);
        probs(ii/interval+1, 4) = numel(find(choice_blur(:, 4) == ii))/size(people, 1);
        probs(ii/interval+1, 5) = numel(find(choice_blur(:, 5) == ii))/size(people, 1);
        probs(ii/interval+1, 6) = numel(find(choice_blur(:, 6) == ii))/size(people, 1);
        probs(ii/interval+1, 7) = sum(probs(ii/interval+1, 1) + probs(ii/interval+1, 2) + ...
            probs(ii/interval+1, 3) + probs(ii/interval+1, 4) +probs(ii/interval+1, 5) + ...
            probs(ii/interval+1, 6));
    end
    figure;
    bar(min_aper:interval:max_aper, probs(:,1:3));
    legend('First choice', 'Second choice', 'Third preference');
    title(names(jj));
    xlabel("Aperture Radii")
    ylabel("Probability")
end
figure;
subplot(2,1,1);
b=bar3(1:SceneNumber, squeeze(scn_p_aper(:, :, 1)));
colormap(pink);
colorbar;
% title("Aperture Preference For Each Person Among All Scenes");
xlabel("Participant Number")
ylabel("Scene Index")
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
std_aperture = std(scn_p_aper(:, :, 1));
mean_aperture = mean(scn_p_aper(:, :, 1));
I = zeros(16, size(people, 1));
for kk=1:size(people, 1)
    if std_aperture(1, kk) == 0
        std_aperture(1, kk) = std_aperture(1, kk) + 0.1;
        I(:, kk) = normpdf((0:15)',mean_aperture(1, kk),std_aperture(1, kk))/5;
    else
        I(:, kk) = normpdf((0:15)',mean_aperture(1, kk),std_aperture(1, kk));
    end
end
[sorted_mean, idx] = sort(mean_aperture);
% [sorted_std, idx] = sort(std_aperture(idx));
I_sorted = I(:, idx);
subplot(2,1,2);
[X,Y] = meshgrid(1:size(people, 1), 0:max_aper);
imagesc(I);
% xticklabels(num2str(idx));
colorbar;
% title("Normal Distribution of Aperture Preference Among Scenes");
xlabel("Participant Number")
ylabel("Aperture")
%%
% compute log-likelihood for three scenes
selected_prob = zeros(SceneNumber, max_aper/interval+1);
log_Lambda = zeros(SceneNumber, max_aper/interval+1, max_aper/interval+1);
for jj=6:SceneNumber
    P = ones(1, size(aperture, 2));
    P = Bradley_Terry_old(W_scenes(:, :, jj), P);
    figure;
    stem(aperture, P, "blue", "filled")
    hold on;
    stem(aperture(4), P(4), "red", "filled")
    hold on;
    if jj ~= 4
        stem(aperture(1), P(1), "red", "filled")
    else
        stem(aperture(2), P(2), "red", "filled")
    end
    hold off;
    xlabel("Aperture Size");
    ylabel("Likelihood");
    xticks(aperture)
    title(names(jj)); 
    xlabel("Aperture Radii")
    ylabel("$q_i$", 'Interpreter', 'latex', 'FontSize', 17)
    ylim([0 5]);
    
    for r1=1:max_aper/interval+1
        for r2=1:max_aper/interval+1
            if r1~=r2
                prob1 = (P(r1))/((P(r1)) + (P(r2)));
                prob2 = 1-prob1;
                n12 = W_scenes(r1, r2, jj);
                n21 = W_scenes(r2, r1, jj);
                % temp = log(0.5^(n12+n21)) - log(prob1^n12 * prob2^n21);
                % if(temp>0)
                    prob1 = n12/(n12+n21);
                    prob2 = 1-prob1;
                    temp = log(0.5^(n12+n21)) - log(prob1^n12 * prob2^n21);
                % end
                log_Lambda(jj, r1, r2) = (-2*temp);
            end
        end
    end
    
    tempLambda = squeeze(log_Lambda(jj, :, :));
    a = sum(log(P));
    for rr=1:6
        selected_prob(jj, rr) = sum(squeeze(W_scenes(rr, :, jj)))...
            /(size(people,1)*comparisonPair);
    end
    
end

%% Bradley-Terry model
function P = Bradley_Terry_Newman(w, P)
    n = size(w, 2);
    idx = 1:n;
    for iter = 1:100
        P_last = P;
        for k = 1:n
            nominator = sum(w(k, idx~=k) .* (P_last(idx~=k) ./ (P_last(k) + P_last(idx~=k))));
            denominator =  sum(w(idx~=k, k)' ./ (P_last(k) + P_last(idx~=k)));
            P(k) = nominator/ denominator;
            
        end
        P = P ./ (prod(P_last)^(1/n));
    end
end
function P = Bradley_Terry_old(w, P)
    n = size(w, 2);
    idx = 1:n;
    for iter = 1:100
        for k = 1:n
            nominator = sum(w(k, idx~=k));
            denominator =  sum((w(idx~=k, k)' + w(k, idx~=k)) ./ (P(k) + P(idx~=k)));
            P(k) = nominator/ denominator;
            if P(k) == 0
                P(k) = 0.01;
            end
            
        end
        P = P ./ (prod(P(P~=0))^(1/n));
    end
end
%% Compute Rankings
function W = computeRankings(people, scene, range, interval)
    W = zeros(6, 6);
    for p = 1:size(people,1)
       strcat("Person", num2str(p))
       person = double(strcmpi(people(p, range), {'B'})) + 1;
       for k = 1:size(scene, 2)
           idx_i = scene(p, k, person(k)) / interval + 1;
           % disp("The Choice")
           % scene(p, k, person(k))
           if person(k) == 1
               idx_j =  scene(p, k, 2) / interval + 1;
               % disp("The Alternative")
               % scene(p, k, 2)
           else
               idx_j =  scene(p, k, 1) / interval + 1;
               % disp("The Alternative")
               % scene(p, k, 1)
           end
            W(idx_i, idx_j) = W(idx_i, idx_j) + 1;
       end
    end
end
function W = computeRankingsSingle(people, scene, range, interval, personID)
   W = zeros(6, 6);
   p = personID;
   person = double(strcmpi(people(p, range), {'B'})) + 1;
   for k = 1:size(scene, 2)
       idx_i = scene(p, k, person(k)) / interval + 1;
       if person(k) == 1
           idx_j =  scene(p, k, 2) / interval + 1;
       else
           idx_j =  scene(p, k, 1) / interval + 1;
       end
        W(idx_i, idx_j) = W(idx_i, idx_j) + 1;
   end
    
end


% Make population vector plots

%% Get population vector correlation with distance - CA1

% tBins = -1.5:dt:1;
% binIDs = (tBins(1:end-1) + tBins(2:end))/2;

colors = cbrewer('qual', 'Set1', 8, 'cubic');
alpha = 0.05;

temp = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 7, 7)), 1);
tempb = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 8, 8)), 1);
% temp = CA1corr(:,:,sub2ind([10 10], 7, 7));
% tempb =  CA1corr(:, :, sub2ind([10 10], 8, 8));
% temp = triu(temp);
% tempb = triu(tempb);
% temp(~(triu(ones(size(temp))))) = nan;
% tempb(~(triu(ones(size(tempb))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% for ii = 2:length(tempb)
%     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% end
temp = [temp; tempb];

temp2 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 9, 9)), 1);
temp2b = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 10, 10)), 1);
% temp2 = CA1corr(:,:,sub2ind([10 10], 9, 9));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end
% temp2b = CA1corr(:,:,sub2ind([10 10], 10, 10));
% temp2b = triu(temp2b);
% temp2b(~(triu(ones(size(temp2b))))) = nan;
% for ii = 2:length(temp2b)
%     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% end
temp2 = [temp2; temp2b];


% test whether they come from the same distributions
p = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([temp(:, ii), temp2(:, ii)], [], 'off');
end

% data = [temp(:); temp2(:)];
% groups = [ones(length(temp(:)), 1); 2*ones(length(temp(:)), 1)];
% [~, p] = kruskalwallis(data, groups);

% Plot correlation distances cond 7 9
f = figure;
f.Position = [10 10 800 800];
subplot(4, 4, [1 ])
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')


box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1.2])
title(['CA1'])

%% CA3 Get population vector correlation with distance
temp = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 7, 7)), 1);
tempb = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 8, 8)), 1);

temp2 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 9, 9)), 1);
temp2b = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 10, 10)), 1);

% temp = CA3corr(:,:,sub2ind([10 10], 7, 7));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% tempb = CA3corr(:,:,sub2ind([10 10], 8, 8));
% tempb = triu(tempb);
% tempb(~(triu(ones(size(tempb))))) = nan;
% for ii = 2:length(tempb)
%     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% end
temp = [temp; tempb];

% temp2 = CA3corr(:,:,sub2ind([10 10], 9, 9));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end
% temp2b = CA3corr(:,:,sub2ind([10 10], 10, 10));
% temp2b = triu(temp2b);
% temp2b(~(triu(ones(size(temp2b))))) = nan;
% for ii = 2:length(temp2b)
%     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% end
temp2 = [temp2; temp2b];

p = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([temp(:, ii), temp2(:, ii)], [], 'off');
end


% data = [temp(:); temp2(:)];
% groups = [ones(length(temp(:)), 1); 2*ones(length(temp(:)), 1)];
% [~, p] = kruskalwallis(data, groups);

% Plot correlation distances cond 1 3 5
figure(f)
subplot(4, 4, [2 ])
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title(['CA3'])
xlim([0 length(binIDs)])
ylim([-0.8 1.2])
%% LS Get population vector correlation with distance
temp = getAvgCorr(LScorr(:,:,sub2ind([10 10], 7, 7)), 1);
tempb = getAvgCorr(LScorr(:,:,sub2ind([10 10], 8, 8)), 1);

temp2 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 9, 9)), 1);
temp2b = getAvgCorr(LScorr(:,:,sub2ind([10 10], 10, 10)), 1);

% temp = LScorr(:,:,sub2ind([10 10], 7, 7));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% tempb = LScorr(:,:,sub2ind([10 10], 8, 8));
% tempb = triu(tempb);
% tempb(~(triu(ones(size(tempb))))) = nan;
% for ii = 2:length(tempb)
%     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% end
temp = [temp; tempb];

% temp2 = LScorr(:,:,sub2ind([10 10], 9, 9));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end
% temp2b = LScorr(:,:,sub2ind([10 10], 10, 10));
% temp2b = triu(temp2b);
% temp2b(~(triu(ones(size(temp2b))))) = nan;
% for ii = 2:length(temp2b)
%     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% end
temp2 = [temp2; temp2b];


p = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([temp(:, ii), temp2(:, ii)], [], 'off');
end

% data = [temp(:); temp2(:)];
% groups = [ones(length(temp(:)), 1); 2*ones(length(temp(:)), 1)];
% [~, p] = kruskalwallis(data, groups);

% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 3)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title(['LS'])
xlim([0 length(binIDs)])
ylim([-0.8 1.2])
%% Compare region 1 corr with other conds

temp = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 9, 1)), 0);
tempb = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 10, 2)), 0);

temp2 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 9, 3)), 0);
temp2b = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 10, 4)), 0);

temp3 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 9, 5)), 0);
temp3b = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 10, 6)), 0);

temp4 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 9, 9)), 1);
temp4b = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 10, 10)), 1);

% 
% temp = CA1corr(:,:,sub2ind([10 10], 9, 1));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% tempb = CA1corr(:,:,sub2ind([10 10], 10, 2));
% tempb = triu(tempb);
% tempb(~(triu(ones(size(tempb))))) = nan;
% for ii = 2:length(tempb)
%     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% end
temp = [temp; tempb];

% temp2 = CA1corr(:,:,sub2ind([10 10], 9, 3));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
%     
% end
% temp2b = CA1corr(:,:,sub2ind([10 10], 10, 4));
% temp2b = triu(temp2b);
% temp2b(~(triu(ones(size(temp2b))))) = nan;
% for ii = 2:length(temp2b)
%     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% end
temp2 = [temp2; temp2b];


% temp3 = CA1corr(:,:,sub2ind([10 10], 9, 5));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end
% temp3b = CA1corr(:,:,sub2ind([10 10], 10, 6));
% temp3b = triu(temp3b);
% temp3b(~(triu(ones(size(temp3b))))) = nan;
% for ii = 2:length(temp3b)
%     temp3b(ii,:) = circshift(temp3b(ii,:), -ii+1);
% end
temp3 = [temp3; temp3b];


% temp4 = CA1corr(:,:,sub2ind([10 10], 9, 9));
% temp4 = triu(temp4);
% temp4(~(triu(ones(size(temp4))))) = nan;
% for ii = 2:length(temp4)
%     temp4(ii,:) = circshift(temp4(ii,:), -ii+1);
% end
% temp4b = CA1corr(:,:,sub2ind([10 10], 10, 10));
% temp4b = triu(temp4b);
% temp4b(~(triu(ones(size(temp4b))))) = nan;
% for ii = 2:length(temp4b)
%     temp4b(ii,:) = circshift(temp4b(ii,:), -ii+1);
% end
temp4 = [temp4; temp4b];

p1 = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p1(ii) = kruskalwallis([temp(:, ii), [temp4(:, ii); nan(size(temp4, 1), 1)]], [], 'off');
end

p = nan(size(temp, 2), 1);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([temp(:, ii), temp2(:, ii), temp3(:, ii)], [], 'off');
end

% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 5)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.alpha      = 0.2;

options.color_line = [0 0 0];
options.color_area = [0.2 0.2 0.2];
options.color_area = options.color_area/max(options.color_area);
plot_areaerrorbar(temp4, options)
hold on
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')
% add significance
inds1 = p1 < alpha;
vals = 1.2*ones(length(p), 1);
hold on
scatter(p_inds(inds1), vals(inds1), 20, [0.3 0.3 0.3], '*')

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA1 9-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1.3])

%% CA3 Compare region 1 corr with other conds
temp = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 9, 1)), 0);
tempb = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 10, 2)), 0);

% temp = CA3corr(:,:,sub2ind([10 10], 9, 1));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% tempb = CA3corr(:,:,sub2ind([10 10], 10, 2));
% tempb = triu(tempb);
% tempb(~(triu(ones(size(tempb))))) = nan;
% for ii = 2:length(tempb)
%     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% end
temp = [temp; tempb];

temp2 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 9, 3)), 0);
temp2b = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 10, 4)), 0);

% temp2 = CA3corr(:,:,sub2ind([10 10], 9, 3));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end
% temp2b = CA3corr(:,:,sub2ind([10 10], 10, 4));
% temp2b = triu(temp2b);
% temp2b(~(triu(ones(size(temp2b))))) = nan;
% for ii = 2:length(temp2b)
%     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% end
temp2 = [temp2; temp2b];

temp3 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 9, 5)), 0);
temp3b = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 10, 6)), 0);

% temp3 = CA3corr(:,:,sub2ind([10 10], 9, 5));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end
% temp3b = CA3corr(:,:,sub2ind([10 10], 10, 6));
% temp3b = triu(temp3b);
% temp3b(~(triu(ones(size(temp3b))))) = nan;
% for ii = 2:length(temp3b)
%     temp3b(ii,:) = circshift(temp3b(ii,:), -ii+1);
% end
temp3 = [temp3; temp3b];

temp4 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 9, 9)), 1);
temp4b = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 10, 10)), 1);

% temp4 = CA3corr(:,:,sub2ind([10 10], 9, 9));
% temp4 = triu(temp4);
% temp4(~(triu(ones(size(temp4))))) = nan;
% for ii = 2:length(temp4)
%     temp4(ii,:) = circshift(temp4(ii,:), -ii+1);
% end
% temp4b = CA3corr(:,:,sub2ind([10 10], 10, 10));
% temp4b = triu(temp4b);
% temp4b(~(triu(ones(size(temp4b))))) = nan;
% for ii = 2:length(temp4b)
%     temp4b(ii,:) = circshift(temp4b(ii,:), -ii+1);
% end
temp4 = [temp4; temp4b];

p1 = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p1(ii) = kruskalwallis([temp(:, ii), [temp4(:, ii); nan(size(temp4, 1), 1)]], [], 'off');
end

p = nan(size(temp, 2), 1);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([temp(:, ii), temp2(:, ii), temp3(:, ii)], [], 'off');
end


% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 6)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.alpha      = 0.2;
options.color_line = [0 0 0];
options.color_area = [0.2 0.2 0.2];
options.color_area = options.color_area/max(options.color_area);
plot_areaerrorbar(temp4, options)
hold on
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);



plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')
% add significance
inds1 = p1 < alpha;
vals = 1.2*ones(length(p), 1);
hold on
scatter(p_inds(inds1), vals(inds1), 20, [0.3 0.3 0.3], '*')
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA3 9-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1.3])
%% LS Compare region 1 corr with other conds

temp = getAvgCorr(LScorr(:,:,sub2ind([10 10], 9, 1)), 0);
tempb = getAvgCorr(LScorr(:,:,sub2ind([10 10], 10, 2)), 0);

temp2 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 9, 3)), 0);
temp2b = getAvgCorr(LScorr(:,:,sub2ind([10 10], 10, 4)), 0);

temp3 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 9, 5)), 0);
temp3b = getAvgCorr(LScorr(:,:,sub2ind([10 10], 10, 6)), 0);

temp4 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 9, 9)), 1);
temp4b = getAvgCorr(LScorr(:,:,sub2ind([10 10], 10, 10)), 1);
% temp = LScorr(:,:,sub2ind([10 10], 9, 1));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% tempb = LScorr(:,:,sub2ind([10 10], 10, 2));
% tempb = triu(tempb);
% tempb(~(triu(ones(size(tempb))))) = nan;
% for ii = 2:length(tempb)
%     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% end
temp = [temp; tempb];

% temp2 = LScorr(:,:,sub2ind([10 10], 9, 3));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end
% temp2b = LScorr(:,:,sub2ind([10 10], 10, 4));
% temp2b = triu(temp2b);
% temp2b(~(triu(ones(size(temp2b))))) = nan;
% for ii = 2:length(temp2b)
%     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% end
temp2 = [temp2; temp2b];

% temp3 = LScorr(:,:,sub2ind([10 10], 9, 5));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end
% temp3b = LScorr(:,:,sub2ind([10 10], 10, 6));
% temp3b = triu(temp3b);
% temp3b(~(triu(ones(size(temp3b))))) = nan;
% for ii = 2:length(temp3b)
%     temp3b(ii,:) = circshift(temp3b(ii,:), -ii+1);
% end
temp3 = [temp3; temp3b];

% temp4 = LScorr(:,:,sub2ind([10 10], 9, 9));
% temp4 = triu(temp4);
% temp4(~(triu(ones(size(temp4))))) = nan;
% for ii = 2:length(temp4)
%     temp4(ii,:) = circshift(temp4(ii,:), -ii+1);
% end
% temp4b = LScorr(:,:,sub2ind([10 10], 10, 10));
% temp4b = triu(temp4b);
% temp4b(~(triu(ones(size(temp4b))))) = nan;
% for ii = 2:length(temp4b)
%     temp4b(ii,:) = circshift(temp4b(ii,:), -ii+1);
% end
temp4 = [temp4; temp4b];

p1 = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p1(ii) = kruskalwallis([temp(:, ii), [temp4(:, ii); nan(size(temp4, 1), 1)]], [], 'off');
end

p = nan(size(temp, 2), 1);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([temp(:, ii), temp2(:, ii), temp3(:, ii)], [], 'off');
end


% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 7)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.alpha      = 0.2;
options.color_line = [0 0 0];
options.color_area = [0.2 0.2 0.2];
options.color_area = options.color_area/max(options.color_area);
plot_areaerrorbar(temp4, options)
hold on
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);


plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')
% add significance
inds1 = p1 < alpha;
vals = 1.2*ones(length(p), 1);
hold on
scatter(p_inds(inds1), vals(inds1), 20, [0.3 0.3 0.3], '*')
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('LS 9-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1.3])
%% CA1 Compare region 1 corr with same dir and opposite dir

temp = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 9, 10)), 0);
temp2 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 7, 8)), 0);
temp3 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 1, 2)), 0);

% temp = CA1corr(:,:,sub2ind([10 10], 9, 10));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end

% temp2 = CA1corr(:,:,sub2ind([10 10], 7, 8));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end

% temp3 = CA1corr(:,:,sub2ind([10 10], 1, 2));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end


% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 9)
hold off
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA1 9-10, 7-8, 1-2')
%% CA3 Compare region 1 corr with same dir and opposite dir
temp = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 9, 10)), 0);
temp2 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 7, 8)), 0);
temp3 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 1, 2)), 0);

% temp = CA3corr(:,:,sub2ind([10 10], 9, 10));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end

% temp2 = CA3corr(:,:,sub2ind([10 10], 7, 8));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end

% temp3 = CA3corr(:,:,sub2ind([10 10], 1, 2));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end


% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 10)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA3 1-1 2 6')

%% LS Compare region 1 corr with same dir and opposite dir

temp = getAvgCorr(LScorr(:,:,sub2ind([10 10], 9, 10)), 0);
temp2 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 7, 8)), 0);
temp3 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 1, 2)), 0);

% temp = LScorr(:,:,sub2ind([10 10], 9, 10));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% 
% temp2 = LScorr(:,:,sub2ind([10 10], 7, 8));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
% end
% 
% temp3 = LScorr(:,:,sub2ind([10 10], 1, 2));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end


% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 11)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('LS 1-1 2 6')

%% Compare region 1 corr with other conds

temp = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 1, 1)), 1);
temp2 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 1, 3)), 0);
temp3 = getAvgCorr(CA1corr(:,:,sub2ind([10 10], 1, 5)), 0);

% 
% temp = CA1corr(:,:,sub2ind([10 10],1, 1));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% % tempb = CA1corr(:,:,sub2ind([10 10], 6, 6));
% % tempb = triu(tempb);
% % tempb(~(triu(ones(size(tempb))))) = nan;
% % for ii = 2:length(tempb)
% %     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% % end
% % temp = [temp; tempb];
% 
% temp2 = CA1corr(:,:,sub2ind([10 10], 1, 3));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
%     
% end
% % temp2b = CA1corr(:,:,sub2ind([10 10], 6, 4));
% % temp2b = triu(temp2b);
% % temp2b(~(triu(ones(size(temp2b))))) = nan;
% % for ii = 2:length(temp2b)
% %     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% % end
% % temp2 = [temp2; temp2b];
% 
% 
% temp3 = CA1corr(:,:,sub2ind([10 10], 1, 5));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end
% temp3b = CA1corr(:,:,sub2ind([10 10], 6, 2));
% temp3b = triu(temp3b);
% temp3b(~(triu(ones(size(temp3b))))) = nan;
% for ii = 2:length(temp3b)
%     temp3b(ii,:) = circshift(temp3b(ii,:), -ii+1);
% end
% temp3 = [temp3; temp3b];



p1 = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p1(ii) = kruskalwallis([[temp(:, ii); nan(size(temp, 1), 1)], temp3(:, ii)], [], 'off');
end

p = nan(size(temp, 2), 1);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([[temp(:, ii); nan(size(temp, 1), 1)], temp2(:, ii)], [], 'off');
end

% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 13)
hold off
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.alpha      = 0.2;

options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
plot_areaerrorbar(temp, options)
hold on
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')
% add significance
inds1 = p1 < alpha;
vals = 1.2*ones(length(p), 1);
hold on
scatter(p_inds(inds1), vals(inds1), 20, [0.3 0.3 0.3], '*')

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA1 1-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1.3])

%% Compare region 1 corr with other conds

temp = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 1, 1)), 1);
temp2 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 1, 3)), 0);
temp3 = getAvgCorr(CA3corr(:,:,sub2ind([10 10], 1, 5)), 0);


% temp = CA3corr(:,:,sub2ind([10 10],1, 1));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% % tempb = CA3corr(:,:,sub2ind([10 10], 6, 6));
% % tempb = triu(tempb);
% % tempb(~(triu(ones(size(tempb))))) = nan;
% % for ii = 2:length(tempb)
% %     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% % end
% % temp = [temp; tempb];
% 
% temp2 = CA3corr(:,:,sub2ind([10 10], 1, 3));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
%     
% end
% % temp2b = CA3corr(:,:,sub2ind([10 10], 6, 4));
% % temp2b = triu(temp2b);
% % temp2b(~(triu(ones(size(temp2b))))) = nan;
% % for ii = 2:length(temp2b)
% %     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% % end
% % temp2 = [temp2; temp2b];
% 
% 
% temp3 = CA3corr(:,:,sub2ind([10 10], 1, 5));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end
% % temp3b = CA3corr(:,:,sub2ind([10 10], 6, 2));
% % temp3b = triu(temp3b);
% % temp3b(~(triu(ones(size(temp3b))))) = nan;
% % for ii = 2:length(temp3b)
% %     temp3b(ii,:) = circshift(temp3b(ii,:), -ii+1);
% % end
% % temp3 = [temp3; temp3b];



p1 = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p1(ii) = kruskalwallis([[temp(:, ii); nan(size(temp, 1), 1)], temp3(:, ii)], [], 'off');
end

p = nan(size(temp, 2), 1);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([[temp(:, ii); nan(size(temp, 1), 1)], temp2(:, ii)], [], 'off');
end

% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 14)
hold off
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.alpha      = 0.2;

options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
plot_areaerrorbar(temp, options)
hold on
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')
% add significance
inds1 = p1 < alpha;
vals = 1.2*ones(length(p), 1);
hold on
scatter(p_inds(inds1), vals(inds1), 20, [0.3 0.3 0.3], '*')

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA3 1-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1.3])

%% Compare region 1 corr with other conds

temp = getAvgCorr(LScorr(:,:,sub2ind([10 10], 1, 1)), 1);
temp2 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 1, 3)), 0);
temp3 = getAvgCorr(LScorr(:,:,sub2ind([10 10], 1, 5)), 0);
% temp = LScorr(:,:,sub2ind([10 10],1, 1));
% temp = triu(temp);
% temp(~(triu(ones(size(temp))))) = nan;
% for ii = 2:length(temp)
%     temp(ii,:) = circshift(temp(ii,:), -ii+1);
% end
% % tempb = LScorr(:,:,sub2ind([10 10], 6, 6));
% % tempb = triu(tempb);
% % tempb(~(triu(ones(size(tempb))))) = nan;
% % for ii = 2:length(tempb)
% %     tempb(ii,:) = circshift(tempb(ii,:), -ii+1);
% % end
% % temp = [temp; tempb];
% 
% temp2 = LScorr(:,:,sub2ind([10 10], 1, 3));
% temp2 = triu(temp2);
% temp2(~(triu(ones(size(temp2))))) = nan;
% for ii = 2:length(temp2)
%     temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
%     
% end
% % temp2b = LScorr(:,:,sub2ind([10 10], 6, 4));
% % temp2b = triu(temp2b);
% % temp2b(~(triu(ones(size(temp2b))))) = nan;
% % for ii = 2:length(temp2b)
% %     temp2b(ii,:) = circshift(temp2b(ii,:), -ii+1);
% % end
% % temp2 = [temp2; temp2b];
% 
% 
% temp3 = LScorr(:,:,sub2ind([10 10], 1, 5));
% temp3 = triu(temp3);
% temp3(~(triu(ones(size(temp3))))) = nan;
% for ii = 2:length(temp3)
%     temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
% end
% % temp3b = LScorr(:,:,sub2ind([10 10], 6, 2));
% % temp3b = triu(temp3b);
% % temp3b(~(triu(ones(size(temp3b))))) = nan;
% % for ii = 2:length(temp3b)
% %     temp3b(ii,:) = circshift(temp3b(ii,:), -ii+1);
% % end
% % temp3 = [temp3; temp3b];



p1 = nan(size(temp, 2), 1);
p_inds = 1:length(p);
for ii = 1:size(temp, 2)
    p1(ii) = kruskalwallis([[temp(:, ii); nan(size(temp, 1), 1)], temp3(:, ii)], [], 'off');
end

p = nan(size(temp, 2), 1);
for ii = 1:size(temp, 2)
    p(ii) = kruskalwallis([[temp(:, ii); nan(size(temp, 1), 1)], temp2(:, ii)], [], 'off');
end

% Plot correlation distances cond 1 3 5
% f = figure;
figure(f)
subplot(4, 4, 15)
hold off
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.alpha      = 0.2;

options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
plot_areaerrorbar(temp, options)
hold on
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)

% add significance
inds = p < alpha;
vals = 1.1*ones(length(p), 1);
hold on
scatter(p_inds(inds), vals(inds), 20, [0.5 0.5 0.5], '*')
% add significance
inds1 = p1 < alpha;
vals = 1.2*ones(length(p), 1);
hold on
scatter(p_inds(inds1), vals(inds1), 20, [0.3 0.3 0.3], '*')

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('LS 1-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1.3])

%% New figure showing the diagonal of control with jump corrs

% Plot correlation diagonal 9 1
f = figure;
f.Position = [10 10 800 800];
subplot(4, 4, [1 ])

plot(diag(CA1corr(:,:,sub2ind([10 10],9, 1))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(CA1corr(:,:,sub2ind([10 10],9, 3))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(CA1corr(:,:,sub2ind([10 10],9, 5))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['CA1'])


%% Plot correlation diagonal 9 1 CA3
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [2 ])

plot(diag(CA3corr(:,:,sub2ind([10 10],9, 1))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(CA3corr(:,:,sub2ind([10 10],9, 3))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(CA3corr(:,:,sub2ind([10 10],9, 5))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['CA3'])


%% Plot correlation diagonal 9 1 LS
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [3 ])

plot(diag(LScorr(:,:,sub2ind([10 10],9, 1))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(LScorr(:,:,sub2ind([10 10],9, 3))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(LScorr(:,:,sub2ind([10 10],9, 5))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['LS'])

%% New figure showing the diagonal of control with jump corrs

% Plot correlation diagonal 10 2
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [5 ])

plot(diag(CA1corr(:,:,sub2ind([10 10],10, 2))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(CA1corr(:,:,sub2ind([10 10],10, 4))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(CA1corr(:,:,sub2ind([10 10],10, 6))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['CA1'])


%% Plot correlation diagonal 9 1 CA3
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [6 ])

plot(diag(CA3corr(:,:,sub2ind([10 10], 10, 2))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(CA3corr(:,:,sub2ind([10 10],10, 4))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(CA3corr(:,:,sub2ind([10 10],10, 6))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['CA3'])


%% Plot correlation diagonal 9 1 LS
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [7 ])

plot(diag(LScorr(:,:,sub2ind([10 10],10, 2))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(LScorr(:,:,sub2ind([10 10],10, 4))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(LScorr(:,:,sub2ind([10 10],10, 6))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['LS'])


%% Average both directions
% Plot correlation diagonal 10 2
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [9 ])
hold off
plot(mean([diag(CA1corr(:,:,sub2ind([10 10],10, 2))) diag(CA1corr(:,:,sub2ind([10 10], 9, 1)))]'), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(mean([diag(CA1corr(:,:,sub2ind([10 10],10, 4))) diag(CA1corr(:,:,sub2ind([10 10], 9, 3)))]'), 'color', colors(2,:), 'lineWidth', 2)
plot(mean([diag(CA1corr(:,:,sub2ind([10 10],10, 6))) diag(CA1corr(:,:,sub2ind([10 10], 9, 5)))]'), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['CA1'])


%% Plot correlation diagonal 9 1 CA3
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [6 ])

plot(diag(CA3corr(:,:,sub2ind([10 10], 10, 2))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(CA3corr(:,:,sub2ind([10 10],10, 4))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(CA3corr(:,:,sub2ind([10 10],10, 6))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['CA3'])


%% Plot correlation diagonal 9 1 LS
figure(f)
f.Position = [10 10 800 800];
subplot(4, 4, [7 ])

plot(diag(LScorr(:,:,sub2ind([10 10],10, 2))), 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(diag(LScorr(:,:,sub2ind([10 10],10, 4))), 'color', colors(2,:), 'lineWidth', 2)
plot(diag(LScorr(:,:,sub2ind([10 10],10, 6))), 'color', colors(3,:), 'lineWidth', 2)


box off
xticks([0:10:40]); xticklabels([0:(dx/10)*10:150])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.3 0.7])
title(['LS'])




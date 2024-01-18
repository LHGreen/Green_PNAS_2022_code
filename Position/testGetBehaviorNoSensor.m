clear all, close all, clc

dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

nchannels = 2;

pos = load([basename '.pos']);

%%
behavior = bz_getJumpBehav_NoPiezoEditVersion(pos);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

% save([basename '.jumpNoPiezo.behavior.mat'])

%%
linear = getBehavEvents_with_Z(pos, 300);

%%

pz = behavior.position.z;
px = behavior.position.x;
py = behavior.position.y;

pxm = movmean(px, 40);
pxrms = fastrms(px, 40);
pxrmsrms = fastrms(fastrms(px, 40), 80);

%%
a = 80000;
b = 82000;
px2 = movmean(px, 40);
figure; plot(px(a:b))
hold on
plot(pxm(a:b))
plot(pxrms(a:b))
plot(pxrmsrms(a:b))
% plot(zscore(px(a:b)))
legend('x', 'movmean(x)', 'rms(x)', 'rms(rms(x))')

%%
vel = diff((px(a:b)));
acc = diff(fastrms(vel, 50));
jerk = diff(acc);
figure
plot(zscore(movmean(px(a:b),10)))
hold on
plot(zscore(movmean(vel, 40)))
plot(zscore(acc))
plot(zscore(jerk))

legend('movmean(x)', 'vel', 'acc', 'jerk')
%%
figure; scatter(px(a:b), py(a:b), 10, 'filled')

%%

difference = behavior.events.jumpTime(:,2) - behavior.events.landTimeNoSensor';
figure; scatter(1:length(difference), difference, 5, 'k', 'filled')
title(['mean = ' num2str(mean(difference)) ', sd = ' num2str(sqrt(var(difference))) ...
    ', CV = ' num2str(sqrt(var(difference))./mean(difference))])
xlabel('trial')
ylabel('sensor time - jerk time (s)')


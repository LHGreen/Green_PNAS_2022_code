
clear all, close all, clc

LFPMat = [];
specMat = [];
specMat2 = [];
phaseMat = [];
condID = [];
sessionID = [];

figDir = fullfile('D:', 'BuzsakiRinzel', 'Figures', 'LFP');
mkdir(figDir);

%%
saveDir = fullfile('D:', 'BuzsakiRinzel', 'Data');

%
thisDir = pwd;

d = dir;

d = d(~ismember({d.name}, {'.', '..'}));

index = zeros(1, length(d));
for ii = 1:length(d)
    index(ii) = isfolder(d(ii).name);
end
d = d(logical(index));

dataDir = fullfile('D:', 'BuzsakiRinzel', 'Data');

region = 'ca1';
type = 'jump';


filterFreq = [6 13]/(1250*0.5); % lfp.samplingRate = 1250
[b, a] = butter(3, filterFreq, 'bandpass');

% try getting spectrogram from each
%%
for ii = 1:length(d)
    
    try
        cd(thisDir)
        cd(d(ii).name)
        fprintf(['Session ' d(ii).name '\n'])
        fprintf([pwd '\n'])
        
        dirName = pwd;
        basename = bz_BasenameFromBasepath(dirName);
        
        load([basename '.behavior.mat']);
        
        behavior.events.jumpTime(1)
        
        sessionInfo = bz_getSessionInfo;
        lfp = bz_GetLFP(sessionInfo.ca1);
        
        phases = angle(hilbert(filtfilt(b, a, double(lfp.data))));
        
        rate = lfp.samplingRate;
        
        %% get good trials
        
        try
            behavior.events.badTrials;
        catch
            badTrials = [];
            % store all LFPs first it's faster
            LFPvals = cell(length(behavior.events.trials), 1);
            for iii = 1:length(behavior.events.trials)
                index = InIntervals(lfp.timestamps, behavior.events.trialIntervals(iii,:));
                LFPvals{iii} = lfp.data(index);
            end
            figure('units', 'normalized', 'outerposition', [0.05 0.05 0.5 0.6])
            for iii = 1:length(behavior.events.trials)
                fprintf(['Trial ' num2str(iii) '/' num2str(length(behavior.events.trials)) '\n'])
                
                plot(LFPvals{iii}, 'k', 'lineWidth', 1)
                box off
                pause
                s = input('Bad trial?', 's');
                if ~isempty(s)
                    badTrials = [badTrials iii];
                end
            end
            clear LFPvals
            close
            behavior.events.badTrials = badTrials;
            dataName = [basename '.behavior.mat'];
            save(dataName, 'behavior')
            
        end
        
        goodTrials = 1:length(behavior.events.trials);
        goodTrials(behavior.events.badTrials) = [];
        
        trials = goodTrials;
        wavespecMat = nan(3125, 150, length(trials));
        wavespecMat2 = nan(3125, 100, length(trials));
        sessionMat = nan(length(trials), lfp.samplingRate*2.5+1);
        phaseSessionMat = nan(length(trials), lfp.samplingRate*2.5+1);
        
        %%
        clear wavespec
        wavespec = bz_WaveSpec(lfp, 'ncyc', 3, 'nfreqs', 150, 'frange', [1 200], ...
            'intervals', [min(behavior.events.trialIntervals(:, 1))-20 behavior.timestamps(end)], 'space', 'lin', 'showprogress', true);
        
        
        LFPtimestamps = linspace(min(behavior.events.trialIntervals(:, 1))-20, behavior.timestamps(end), length(wavespec.data));
        
        
        for iii = 1:length(trials)
            if behavior.events.trialConditions(trials(iii)) < 7
                time = behavior.events.jumpTime(trials(iii),1);
            else
                time = round(mean(behavior.events.trialIntervals(trials(ii),:)), 3);
            end
            index = InIntervals(LFPtimestamps, [time-1.5 time+1]);
            temp = wavespec.data(index,:);
            if length(temp)~=3125
                ind1 = find(index == 1, 1);
                temp = wavespec.data(ind1:ind1 + 3124,:);
            end
            
            wavespecMat(:,:,iii) = temp;
        end
        
        wavespec1Freqs = wavespec.freqs;
        
        %%
        clear wavespec
        wavespec = bz_WaveSpec(lfp, 'ncyc', 3, 'nfreqs', 100, 'frange', [1 50], ...
            'intervals', [min(behavior.events.trialIntervals(:, 1))-20 behavior.timestamps(end)], 'space', 'lin', 'showprogress', true);
        
        for iii = 1:length(trials)
            if behavior.events.trialConditions(trials(iii)) < 7
                time = behavior.events.jumpTime(trials(iii),1);
            else
                time = round(mean(behavior.events.trialIntervals(trials(ii),:)), 3);
            end
            index = InIntervals(LFPtimestamps, [time-1.5 time+1]);
            temp = wavespec.data(index,:);
            if length(temp)~=3125
                ind1 = find(index == 1, 1);
                temp = wavespec.data(ind1:ind1 + 3124,:);
            end
           
            wavespecMat2(:,:,iii) = temp;
        end
        
        
        for iii = 1:length(trials)%1:length(behavior.events.trials)
            ntrial = trials(iii);
            if behavior.events.trialConditions(ntrial) < 7
                index = round(rate*behavior.events.jumpTime(ntrial, 2));
            else
                index = round(rate*mean(behavior.events.trialIntervals(ntrial, :)));
            end
            if ~isnan(index)
                sessionMat(iii, :) = lfp.data((index-rate*1.5):(index+rate*1));
                phaseSessionMat(iii, :) = phases((index-rate*1.5):(index+rate*1));
            end
        end
        %%
        
        LFPMat = [LFPMat; sessionMat];
        specMat = cat(3, specMat, wavespecMat);
        specMat2 = cat(3, specMat2, wavespecMat2);
        phaseMat = [phaseMat; phaseSessionMat];
        condID = [condID; behavior.events.trialConditions(trials)'];
        sessionID = [sessionID; ii*ones(length(trials), 1)];
        
        fprintf([basename ' added \n'])
       
        
    catch
        
        fprintf([basename ' did not work \n'])
    end
    
    cd ..
    
    clear lfp
end

colors = cbrewer('qual', 'Set2', 8, 'cubic');
%%

dataName = ['LFP_Spec' type 'FINAL_2' region '.mat'];

cd('LFPData')
save(dataName, 'specMat', 'specMat2', 'LFPMat', 'phaseMat', 'condID', 'sessionID', 'd', '-v7.3')

cd(thisDir)

%% test specMat
jump = condID < 7;
run = condID > 6;


medianWavespec = median((abs(specMat(:,:,jump))), 3);
medianWavespec = medianWavespec.^2;

medianWavespec_run = median((abs(specMat(:,:,run))), 3);
medianWavespec_run = medianWavespec_run.^2;
%%
fig = figure;
subplot(5, 2, 1)
set(gcf, 'Position', [50 50 500 800])
imagesc(log(medianWavespec'))
axis xy
yticks(38:38:150); yticklabels(50:50:175)
xticks(1:625:3125); xticklabels(-1.5:0.5:1)
colorbar
% clim(1*[-1 1])
box off
xlabel('time (s)')
ylabel('Frequency (Hz)')
title(['nTrials = ' num2str(sum(jump)) ])

subplot(5, 2, 3)
imagesc(log(medianWavespec_run'))
axis xy
yticks(38:38:150); yticklabels(50:50:175)
xticks(1:625:3125); xticklabels(-1.5:0.5:1)
colorbar
% clim(1*[-1 1])
box off
xlabel('time (s)')
ylabel('Frequency (Hz)')
title(['nTrials = ' num2str(sum(run)) ])
%%
figname = [region type '_logLFP.png'];
oldDir = cd(figDir);
saveas(gcf, figname)
cd(oldDir)
close

%% Plot for all conditions
fig = figure;
fig.Position = [10 10 500 800];
for cond = 1:8
    inds = condID == cond;
    
    medianWavespec = median(abs(specMat(:, :, inds)), 3);
    medianWavespec = medianWavespec.^2;
    
    subplot(5, 2, cond)
    imagesc(log(medianWavespec'))
    axis xy
    yticks(38:38:150); yticklabels(50:50:175)
    xticks(1:625:3125); xticklabels(-1.5:0.5:1)
    colorbar
    box off
    %     clim([8 17])
    xlabel('time (s)')
    ylabel('Frequency (Hz)')
    title(num2str(cond))
    
end


%% test specMat
medianWavespec2 = median((abs(specMat2(:,:,jump))), 3);
medianWavespec2 = medianWavespec2.^2;
medianWavespec2_run = median((abs(specMat2(:,:,run))), 3);
medianWavespec2_run = medianWavespec2_run.^2;
%% normalize by median power (similar effect to z-scoring though)
medianWavespec2 = (medianWavespec2./median(medianWavespec2));

%% normalize by baseline

medianWavespec2 = medianWavespec2./median(medianWavespec2(end:-1:end-10, :));
%%

fig = figure;
subplot(5, 2, 1)
set(gcf, 'Position', [50 50 500 800])
imagesc(log(medianWavespec2'))
axis xy
yticks(20:20:100); yticklabels(10:10:50)
xticks(1:625:3125); xticklabels(-1.5:0.5:1)
% clim([13.5 19])
box off
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title(['nTrials = ' num2str(sum(jump)) ])

subplot(5, 2, 3)
imagesc(log(medianWavespec2_run'))
axis xy
yticks(20:20:100); yticklabels(10:10:50)
xticks(1:625:3125); xticklabels(-1.5:0.5:1)
% clim([13.5 17.5])
box off
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title(['nTrials = ' num2str(sum(run)) ])
%%
colors = cbrewer('qual', 'Set2', 8);
fig = figure;
fig.Position = [10 10 500 800];
subplot(5, 2, 1)

% 1.5 seconds before jump
plot(log(medianWavespec2(1, 1:50)), 'lineWidth', 2)
hold on
% at jump
plot(log(medianWavespec2(1875, 1:50)), 'lineWidth', 2)
% 200 ms after jump
plot(log(medianWavespec2(2250,1:50)), 'lineWidth', 2)
xticks(20:20:50); xticklabels(10:10:25)
xlim([1 50])
ylabel('log(power)')
xlabel('Frequency')
box off

%% Okay so for each session, get median spectrogram and find the peaks
sessions = unique(sessionID);
sessions(3) = [];
pre = nan(length(sessions), 100);
takeoff = nan(length(sessions), 100);
land = nan(length(sessions), 100);
control = nan(length(sessions), 100);


for ii = 1:length(sessions)
    pre(ii, :) = median(abs(specMat2(5, :, jump & sessionID == sessions(ii))),3).^2;
    takeoff(ii, :) = median(abs(specMat2(1875, :, jump & sessionID == sessions(ii))),3).^2;
    land(ii, :) = median(abs(specMat2(2250, :, jump & sessionID == sessions(ii))),3).^2;
    control(ii, :) = median(abs(specMat2(1875, :, run & sessionID == sessions(ii))),3).^2;
end

%%
colors = cbrewer('qual', 'Set2', 8);
fig = figure;
fig.Position = [10 10 500 800];
subplot(5, 2, 1)
options.handle = fig;
options.error = 'sem';
options.line_width = 2;
options.alpha = 0.2;
options.color_line = colors(1,:);
options.color_area = colors(1,:)*1.5;
options.color_area = options.color_area/max(options.color_area);
% 1.5 seconds before jump
plot_areaerrorbar(log(pre), options)
hold on
% at jump
options.color_line = colors(2,:);
options.color_area = colors(2,:)*1.5;
options.color_area = options.color_area/max(options.color_area);
% 1.5 seconds before jump
plot_areaerrorbar(log(takeoff), options)
hold on
% 200 ms after jump
options.color_line = colors(3,:);
options.color_area = colors(3,:)*1.5;
options.color_area = options.color_area/max(options.color_area);
% 1.5 seconds before jump
plot_areaerrorbar(log(land), options)
hold on
options.color_line = colors(4,:);
options.color_area = colors(4,:)*1.5;
options.color_area = options.color_area/max(options.color_area);
% 1.5 seconds before jump
plot_areaerrorbar(log(control), options)
xticks(20:20:50); xticklabels(10:10:25)
xlim([1 50])
ylabel('log(power)')
xlabel('Frequency')
box off

%% Then get max of each to plot boxplots of theta frequency

[~, preInd] = max(log(pre(:, 10:40))');
[~, takeoffInd] = max(log(takeoff(:, 10:40))');
[~, landInd] = max(log(land(:, 10:40))');
[~, controlInd] = max(log(control(:, 10:40))');

freqs = linspace(1, 50, 100);

data = nan(length(sessions), 4);
data(:, 1) = freqs(preInd + 10);
data(:, 2) = freqs(takeoffInd + 10);
data(:, 3) = freqs(landInd + 10);
data(:, 4) = freqs(controlInd + 10);

[p, ~, stats] = kruskalwallis(data);
multcompare(stats)

fig = figure;
fig.Position = [10 10 500 800];
subplot(5, 2, 1)
boxplot(data, 'PlotStyle', 'compact')
box off
ylabel('Theta Frequency (Hz)')
ylim([5 15])% ylim([5 13])
subplot(5, 2, 5)
title(num2str(mean(data)))
subplot(5, 2, 7)
title(num2str(std(data)))
%%
figname = [region type '_logLFP2.png'];
oldDir = cd(figDir);
saveas(gcf, figname)
cd(oldDir)
close


%% Now get the phases and plot the probability of phase

filterFreq = [6 12]/(lfp.samplingRate*0.5);
[b, a] = butter(3, filterFreq, 'bandpass');
phases = nan(size(LFPMat));
for ii = 1:size(LFPMat, 1)
    phases(ii,:) = angle(hilbert(filtfilt(b, a, double(LFPMat(ii,:)))));
end


%% Plot the phase probability heat map
dphi = 0.1;
phaseBins = -pi:dphi:pi;
nPhiBins = length(phaseBins)-1;

dt = 0.015;
timeBins = -1.5:dt:1;

sd = 0.7;

time = linspace(-1.5, 1, 3126);

jt = find(condID < 7);

phaseJump = phaseMat(jt, :);

p_phase = histcounts2(reshape(phaseJump', 1, []), repmat(time, 1, size(phaseJump, 1)), phaseBins, timeBins, 'Normalization', 'probability');

fig = figure;
subplot(5, 2, 1)
imagesc(imgaussfilt(p_phase, sd))
ylabel('phase')
xlabel('time')
xticks(1:0.5/dt:length(timeBins)-1); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
clim([0 1.5e-4])
box off
title('Jump')

jt = find(condID > 6);

phaseJump = phaseMat(jt, :);

p_phase = histcounts2(reshape(phaseJump', 1, []), repmat(time, 1, size(phaseJump, 1)), phaseBins, timeBins, 'Normalization', 'probability');


subplot(5, 2, 3)
imagesc(imgaussfilt(p_phase, sd))
ylabel('phase')
xlabel('time')
xticks(1:0.5/dt:length(timeBins)-1); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
clim([0 1.5e-4])
box off
title('Run')
fig.Position = [20 20 1000 750];

%%
figname = [region type '_phaseProb.png'];

oldDir = cd(figDir);
saveas(gcf, figname)
cd(oldDir)
close

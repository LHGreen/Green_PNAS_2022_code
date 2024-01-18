% Go through behavior and refine landing time
clear all, close all; clc
dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

nchannels = 2;
% save here
saveDir = fullfile('D:', 'BuzsakiRinzel', 'Data', basename);

try
load([basename '.behavior.mat'])
catch
    cd(saveDir)
    load([basename '.behavior.mat'])
    cd(dirName)
end

% make new vectors for land time and position
behavior.events.jumpTimePiezo = behavior.events.jumpTime;
behavior.events.jumpLocPiezo = behavior.events.jumpLoc;

landLocs = behavior.events.jumpTime(:, 2);

% for each landing location, find 2 seconds before, 1 second after, and get
% sensor landing time
startTime = behavior.timestamps(1);
startInd = round(behavior.timestamps(1)*120);
%%

figure('units', 'normalized', 'outerposition', [0 0 1 1])
ii = 0;
for time = landLocs'
    ii = ii + 1;
    
    if ~isnan(time)
        sensor = bz_LoadBinary('analogin.dat', 'nchannels', nchannels, ...
            'start', time-1, 'duration', 2);
        
        
        sensor = (movmean(sensor, 40));
        [max1, id1] = max(abs(zscore(sensor(:,1))));
        [max2, id2] = max(abs(zscore(sensor(:,nchannels))));
        
        % The first sensor is probably on the left side. Try this.
        
        if mod(behavior.events.trialConditions(ii), 2) == 0
            landSensor = 1;
        else
            landSensor = 2;
        end
        
        meanSensor = mean(sensor(1:10e3,landSensor));
        sdSensor = sqrt(var(sensor(1:10e3,landSensor)));
        sensorLandTs = find(abs(sensor(1:end,landSensor)-meanSensor) > 6*sdSensor,1); % get landing time
%         sensorLandTs = sensorLandTs + 20e3;
        % convert to overall time index
        % convert to 120 Hz, add previous landing time, and
        % subtract the 2 s buffer
        landTs = round(sensorLandTs/167+time*120 - 120)-startInd;
        % jump Index
        jumpTs = round(behavior.events.jumpTime(ii, 1)*120)-startInd;
        
        timeIndex = round(time*120)-startInd;
        % Then double check landing time
        
        % plot
        
        
        subplot(2,2,1)
        if mod(behavior.events.trialConditions(ii), 2) == 0
            title('even')
        else
            title('odd')
        end
        plot(behavior.position.x(timeIndex -120: timeIndex + 120), behavior.position.z(timeIndex -120: timeIndex + 120))
        hold on
        scatter(behavior.position.x(jumpTs), behavior.position.z(jumpTs), 30, 'g', 'filled')
        scatter(behavior.position.x(landTs), behavior.position.z(landTs), 30, 'r', 'filled')
        scatter(behavior.position.x(timeIndex), behavior.position.z(timeIndex), 30, 'b', 'filled')
        legend('xz position', 'jump', 'sensor land', 'acc land', 'location', 'southeast')
        xlabel('x')
        ylabel('z')
        box off
        if mod(behavior.events.trialConditions(ii), 2) == 0
            title('even')
        else
            title('odd')
        end
        hold off
        
        subplot(2, 2, 2)
        plot(0, 0)
        title(['trial ' num2str(ii) '/' num2str(sum(~isnan(landLocs)))])
        
        subplot(2,2,4)
        plot(behavior.position.x(timeIndex -120: timeIndex + 120))
        hold on
        scatter(jumpTs-timeIndex+120, behavior.position.x(jumpTs), 30, 'g', 'filled')
        scatter(landTs-timeIndex + 120, behavior.position.x(landTs), 30, 'r', 'filled')
        scatter(121, behavior.position.x(timeIndex), 30, 'b', 'filled')
        legend('x position', 'jump', 'sensor land', 'acc land', 'location', 'northeast')
        xlabel('time')
        ylabel('x')
        box off
        hold off
        
        subplot(2,2,3)
        plot(zscore(sensor(:,1)))
        hold on
        plot(zscore(sensor(:,2)))
        scatter(-diff(behavior.events.jumpTime(ii,:))*20e3+20e3, 0, 30, 'g', 'filled')
        scatter(sensorLandTs, zscore(sensor(sensorLandTs, landSensor)), 30, 'r', 'filled')
        legend('s1', 's2', 'jump', 'land', 'location', 'northeast')
        sensorLandTs
        xlabel('time')
        ylabel('sensor value')
        box off
        hold off
        shg
        pause; %s{k}(i) = 'y';
        prompt = input('Throw out: [Y]','s');
        
        % If good
        if isempty(prompt)
            landTime = behavior.timestamps(landTs);
        else
            landTime = nan;
            landTs = nan;
        end
    end
    
    
    % take care of linear trials
    if isnan(time)
        landTime = nan;
        landTs = nan;
    end
    
    % Then add landTs in behavior file and position
    behavior.events.jumpTimePiezo(ii,2) = landTime;
    if isnan(landTs)
        behavior.events.jumpLocPiezo(ii,:) = nan;
    else
        behavior.events.jumpLocPiezo(ii,[2 4 6]) = [behavior.position.x(landTs) behavior.position.y(landTs) behavior.position.z(landTs)];
    end
    
end
close

%% Show changes
figure; 
subplot(2,2,1)
plot(diff(behavior.events.jumpTime'))
hold on; plot(diff(behavior.events.jumpTimePiezo'))
ylabel('time')
xlabel('trial num')
box off

subplot(2,2,2)
plot(diff(behavior.events.jumpTime')-diff(behavior.events.jumpTimePiezo'))
ylabel('time difference')
xlabel('trial num')
box off

subplot(2, 2, 3) 
plot(behavior.events.jumpLoc(:, 2))
hold on; plot(behavior.events.jumpLocPiezo(:, 2))
ylabel('x position')
xlabel('trial num')
box off

%% Then go save it
cd(saveDir)
filename = [basename '.behavior.mat'];
save(filename, 'behavior')
cd(dirName)

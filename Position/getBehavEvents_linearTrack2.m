function [behavior] = getBehavEvents_linearTrack2(pos,trackEnds)
dbstop if error

if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;
end

% pos(:,8:10) = pos(:,8:10)*1000;

behavior = pos2behav(pos, 'behavType', 'linearTrack');

% idx = intersect(find(pos(:,8) < 2000), find(pos(:,8) > -1400));
scatter(pos(:,8),pos(:,10),'.')
% xlim([0.7 1.7])
% ylim([0.5 2.5])
% ylim([min(pos(:,10)) max(pos(:,10))])
% xlim([min(pos(:,8)) max(pos(:,8))])
ylim([-1200 1000])
xlim([-700 1800])
% axis([-1 1 -1 1])
title('label ends of tracks to get trials')
[startx starty]=ginput();

cx =  startx(end-1:end);
starty = starty(end-1:end);
% Plot distance to ends of track
dist(:,1)=fastrms(abs(pos(:,8)-startx(1))+abs(pos(:,10)-starty(1)),120);
dist(:,2)=fastrms(abs(pos(:,8)-startx(2))+abs(pos(:,10)-starty(2)),120);

plot(dist,'.')
% ylim([-0.2 4])
ylim([-50 2500])

% axis([0 length(dist) 0 3])
[xx yy] = ginput();
close

% So this finds the ends of the track? And choose a minimum distance for it
% to be a real trial? So like about half the track?
[pks_one locs_one]=findpeaks(-dist(:,1),'MINPEAKHEIGHT',-yy, 'MinPeakDistance', 80);
[pks_two locs_two]=findpeaks(-dist(:,2),'MINPEAKHEIGHT',-yy, 'MinPeakDistance', 80);

c=1;
lastTrial = 0;
for i=1:length(locs_one) % for each trial end
    f = locs_one(i)-locs_two; % difference in timestep between ends of trials
    ff = find(f>0); % get all ends before
    [a b]=min(locs_one(i)-locs_two(ff)); % Get first end before, a: distance between ends; b: index of second end
    if a < 1300 & a > 60 & locs_two(ff(b)) > lastTrial + 40 % if less than ~4 s and greater than 0.5s between ends, and greater than last trial by 0.33s
        trials{c} = pos(locs_two(ff(b)):locs_one(i),:);
        c=1+c;
        lastTrial = locs_one(i); %locs_two(ff(b));
    end
    f = locs_two-locs_one(i); 
    ff = find(f>0); % get all ends after
    [a b]=min(locs_two(ff)-locs_one(i)); % get first end after
    if a < 1300 & a > 60 & locs_one(i) > lastTrial + 40 % same criteria as above
        trials{c} = pos(locs_one(i):locs_two(ff(b)),:);
        c=1+c;
        lastTrial = locs_two(ff(b));
    end
    
end

%% merge trials to the same length

trials_unsorted = trials;


size(trials_unsorted)
pause

%%
startx = trackEnds(:,1);
starty = trackEnds(:,2);
endx = trackEnds(:,3);
endy = trackEnds(:,4);
trials = cell(2,1);
b = 1;
% sort trials by direction and shorten to actual trial lengths
for ii = 1:length(trials_unsorted)
    % get x and y position
    xPos = movmean(trials_unsorted{ii}(:,8), 5);
    yPos = movmean(trials_unsorted{ii}(:,10), 5);
    
    % if in direction 1
    if sum(diff(trials_unsorted{ii}(:, 8))) > 0
        % gets distance between start points and all other points in trial
        startDist = pdist2([startx(1) starty(1)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
        endDist = pdist2([endx(2) endy(2)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
        % find index of closest point
        [pks1, start_point] = findpeaks(-startDist, 'MinPeakHeight', -90, 'MinPeakWidth', 8);
        [pks2, end_point] = findpeaks(-endDist, 'MinPeakHeight', -90, 'MinPeakWidth', 8);
        
        if ~isempty(start_point) && ~isempty(end_point)
            start_point = start_point(end)+b;
            end_point = end_point(1) +b;
            
            trials{1}{end+1} = trials_unsorted{ii}(start_point:end_point, :);
        end
    elseif sum(diff(trials_unsorted{ii}(:,8))) < 0
        startDist = pdist2([startx(2) starty(2)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
        endDist = pdist2([endx(1) endy(1)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
        [pks1, start_point] = findpeaks(-startDist, 'MinPeakHeight', -90, 'MinPeakWidth', 1);
        [pks2, end_point] = findpeaks(-endDist, 'MinPeakHeight', -90, 'MinPeakWidth', 1);
        
        if ~isempty(start_point) && ~isempty(end_point)
            start_point = start_point(end)+b;
            end_point = end_point(1) +b;
            
            trials{2}{end+1} = trials_unsorted{ii}(start_point:end_point, :);
        end
    end
end
% Okay I guess that's it?

% Run the checks here to make sure jump trials aren't sorted in
% Then put into behavior struct

% remove trials with low Z values
for condition = 1:length(trials)
    zMin = [];
    for ii = 1:length(trials{condition})
        zMin(ii) = min(trials{condition}{ii}(30:end-30, 9));
    end
    edges = 600:5:750;
    figure; histogram(zMin, edges); shg
    threshold = input('What is the z - threshold between jump and non-jump trials?'); close
    notLinear = find(zMin < threshold);
    for ii =notLinear
        trials{condition}{ii} = [];
    end
    trials{condition}(cellfun('isempty', trials{condition})) = [];
end

% remove all trials with z-max greater than 800
for condition = 1:length(trials)
    for ntrial = 1:length(trials{condition})
        if max(trials{condition}{ntrial}(:,9)) > 800
            trials{condition}{ntrial} = [];
        end
    end
    trials{condition}(cellfun('isempty', trials{condition})) = [];
end

trialCount = 0;
for cond = 1:length(trials)
    for ntrial = 1:length(trials{cond})
        trialCount = trialCount + 1;
        behavior.events.trials{trialCount}.x = trials{cond}{ntrial}(:,8);
        behavior.events.trials{trialCount}.y = trials{cond}{ntrial}(:,10);
        behavior.events.trials{trialCount}.z = trials{cond}{ntrial}(:,9);
        behavior.events.trials{trialCount}.mapping = nan;
        behavior.events.trials{trialCount}.timestamps = trials{cond}{ntrial}(:,1);
        behavior.events.trials{trialCount}.errorPerMarker = trials{cond}{ntrial}(:,11);
        behavior.events.trials{trialCount}.orientation.x = trials{cond}{ntrial}(:,4);
        behavior.events.trials{trialCount}.orientation.y = trials{cond}{ntrial}(:,5);
        behavior.events.trials{trialCount}.orientation.z = trials{cond}{ntrial}(:,6);
        behavior.events.trials{trialCount}.orientation.w = trials{cond}{ntrial}(:,7);
        
        behavior.events.trialConditions(trialCount) = cond;
        behavior.events.trialIntervals(trialCount,:) = [trials{cond}{ntrial}(1,1) trials{cond}{ntrial}(end, 1)];
        
    end
end

behavior.events.conditionType = 'linearTrack';
behavior.events.startEndPos = [startx starty endx endy];



end
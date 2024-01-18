function [behavior] = bz_getJumpBehav_piezoEditVersion(pos, nchannels)

% this does not work well.

dbstop if error

if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;
end
p = pos(:,[8 10]);
behavior = pos2behav(pos,'behavType','jump');

p(:,1) = movmean(p(:,1),40);
p(:,2) = movmean(p(:,2),40);

[r1, r2, r3] = quat2angle([pos(:,7), pos(:,5), pos(:,4), pos(:,6)]);

r2 = circ_smoothTS(r2, 5, 'method', 'mean');
r = diff(r2);

figure
% scatter(pos(:,8),pos(:,9),'.')
scatter(pos(:,8)-pos(:,10),pos(:,9),'.') % x-z is horizontal direction, y is height
% axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
axis([-300 900 400 1000])%([-1000 2000 0 1000])
[x, y] = ginput(); % get click input
close

% Then get start/stop positions
figure
scatter(pos(:,8),pos(:,10),'.') % x-z is horizontal direction, y is height
title('start points')
[startx, starty] = ginput(2);
close

figure
scatter(pos(:,8),pos(:,10),'.') % x-z is horizontal direction, y is height
[endx, endy] = ginput(2);
title('end points')
close

for i=1:length(x) % for all mouse clicks
    dists_x(i,:) = abs(pos(:,8)-pos(:,10)-x(i));
    dists_y(i,:) = abs(pos(:,9)-y(i));
    
    % find points in pos that are closest to mouse click
    [pks_x{i} locs_x{i}] = findpeaks(-dists_x(i,:),'MinPeakHeight',-500);
    [pks_y{i} locs_y{i}] = findpeaks(-dists_y(i,:),'MinPeakHeight',-50);
    for j=1:length(x)
        trials{i,j,1}=[];
        trials{i,j,2}=[];
    end
    jump_ts{i} = 0;
    % for every peak that is close to a mouse click
    for j=length(pks_x{i}):-1:1
        for k=length(pks_y{i}):-1:1
            if abs(pks_x{i}(j)) < 150 && abs(pks_y{i}(k)) < 20 % check that point is reasonably close
                if abs(locs_x{i}(j) - locs_y{i}(k)) < 10 % if aligned in time
                    if abs(mean([locs_x{i}(j) locs_y{i}(k)]) - jump_ts{i}(end)) > 150 % if at least 150 units after previous jump_ts
                        jump_ts{i} = [jump_ts{i};round(mean([locs_x{i}(j) locs_y{i}(k)]))]; % jump ts is mean of x and y indices
                    end
                end
            end
        end
    end
    jump_ts{i} = jump_ts{i}(2:end);
end


count = 1;
% loop through all potential jump times
for k = 1:length(x)
    for i = 1:length(jump_ts{k})
        if jump_ts{k}(i) > 1
            
            time = pos(jump_ts{k}(i),1);
            
            % get piezo sensor info for the trial
            sensor = bz_LoadBinary('analogin.dat', 'nchannels', nchannels, ...
                'start', time-1, 'duration', 5);
            
            sensor = (movmean(sensor, 40));
            
            % Find jump sensor by max sensor value
            [~, id1] = max(zscore(sensor(:,1)));
            [~, id2] = max(zscore(sensor(:,2)));
            
            %             sensorJumpTs = 0;
            %             sensorMax = 0;
            
            landSensor = 0;
            if id1 < id2
                landSensor = 2;
            elseif id1 > id2
                landSensor = 1;
            end
            
            if landSensor ~=0
                if jump_ts{k}(i) < size(pos, 1) -1000
                    meanSensor = mean(sensor(8e3:10e3,landSensor));
                    sdSensor = sqrt(var(sensor(8e3:10e3,landSensor)));
                    sensorLandTs = find(abs(sensor(:,landSensor)-meanSensor) > 7*sdSensor,1);
                    landTs = (round((sensorLandTs-1*20e3)/167)+jump_ts{k}(i));
                    % this is landTs in terms of behavioral timestamps
                    
                    % landTs relative to px
                    landTsTrial = landTs-jump_ts{k}(i)+120;
                    
                    fprintf(['jump_ts{k}{i}] = ' num2str(jump_ts{k}(i)) '\n'])
                    p = round(pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,9), 3); % z-position
                    px = round(pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,8),3); % x-position
                    % x-position variables
                    vel = smoothdata(diff(smoothdata(px, 'sgolay', 5)),'sgolay', 5);
                    acc = smoothdata(diff(vel), 'sgolay', 5);
                    jerk = diff(diff(fastrms(fastrms(diff(movmean(px,40)), 40),80)));
                    jerk2 = diff(acc);
                    
                    rot = r(jump_ts{k}(i)-120:jump_ts{k}(i)+480);
                    rval = r2(jump_ts{k}(i)-120:jump_ts{k}(i)+480);
                    
                    % z-position variables
                    pz = smoothdata(p, 'gaussian', 10);
                    velz = smoothdata(diff(pz),'gaussian', 10);
                    accz = movmean(diff(velz), 10);
                    jerkz = diff(accz);
                    
                    %                     meanPx = mean(px(50:100));
                    %                     stdPx = std(px(50:100));
                    %                     jumpTs = find(abs(px(100:end-100)-meanPx) > 3*stdPx);
                    %
                    %                     fprintf(['jumpTs = ' num2str(jumpTs(1))])
                    %                     jumpTs = jumpTs(1) + 100;
                    
                    if ~isempty(landTsTrial)
                        fprintf(['landTsTrial = ' num2str(landTsTrial)])
                        
                        if landTsTrial < 60
                            [~, jumpTs] = max((rval(1:landTsTrial)));
                        else
                            [~, jumpTs] = max((rval((landTsTrial-60):landTsTrial)));
                            %
                            jumpTs = jumpTs + landTsTrial-61;
                        end
                        
                        
                        if ~isempty(jumpTs)
                            %                     if (jumpTs+jump_ts{k}(i)-120) < landTs % if jump time is before land time
                            % if jump time is at least 0.1 ms and no greater
                            % than 0.5 ms
                            %                         if abs((jumpTs+jump_ts{k}(i)-120)-landTs) < 200 && abs((jumpTs+jump_ts{k}(i)-120)-landTs) > 12
                            
                            % Start by plotting the putative points
                            % plot position in horizontal plane
                            subplot(2,2,2)
                            %                             scatter(pos(:,8),pos(:,10),'.')
                            %                             hold on
                            %                             % plot position around jump time only
                            %                             scatter(pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,8),pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,10),'.r')
                            %                             title(['count = ', num2str(count)])
                            %                             hold off
                            plot(linspace(1,10e4, length(pz)), zscore(pz))
                            hold on
                            plot(linspace(1,10e4, length(velz)), zscore(velz))
                            plot(linspace(1,10e4, length(accz)), zscore(accz))
                            plot(linspace(1,10e4, length(jerkz)), zscore(jerkz))
                            scatter(jumpTs*167, zscore(jerk(jumpTs)), 20, 'g', 'filled')
                            %back to normal
                            scatter(sensorLandTs, zscore(sensor(sensorLandTs, landSensor)), 20, 'r', 'filled')
                            title('z')
                            hold off
                            
                            subplot(2,2,1)
                            % plot position in the xy direction
                            %                             scatter(pos(:,8),pos(:,9),'.')
                            %                             hold on
                            %                             % plot position around the jump time only
                            %                             scatter(pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,8),pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,9),'.r')
                            %                             axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
                            %                             hold off
%                             plot(linspace(1,10e4, length(jerk)), zscore(jerk))
                            
                            scatter(jumpTs*167, zscore(jerk(jumpTs)), 20, 'g', 'filled')
                            hold on
                            scatter(sensorLandTs, zscore(sensor(sensorLandTs, landSensor)), 20, 'r', 'filled')
                            plot(linspace(1,10e4, length(px)), px/max(abs(px)))
                            plot(linspace(1, 10e4, length(rot)), rot/max(rot))
                            plot(linspace(1, 10e4, length(rot)), circ_smoothTS(r2(jump_ts{k}(i)-120:jump_ts{k}(i)+480), 5, 'method', 'mean'))
                            plot(linspace(1, 10e4, length(rot)), r2(jump_ts{k}(i)-120:jump_ts{k}(i)+480))
                            plot(linspace(1, 10e4, length(rot)), (r2(jump_ts{k}(i)-120:jump_ts{k}(i)+480)+ pz)/600)
                            plot(linspace(1, 10e4, length(rot)), zeros(1, length(rot)), 'k')
%                             plot(linspace(1,10e4, length(vel)), vel/max(vel))
%                             plot(linspace(1,10e4, length(acc)), acc/max(acc))
%                             plot(linspace(1, 10e4, length(jerk2)), jerk2/max(jerk2))
                            legend('jumpTime', 'landTime', 'px', 'rot')
                            hold off
                            title('x')
                            
                            
                            % plot sensor input
                            subplot(2,2,3)
                            plot(zscore(sensor(:,1)))
                            hold on
                            plot(zscore(sensor(:,2)))
                            scatter(jumpTs*167, zscore(jerk(jumpTs)), 20, 'g', 'filled')
                            scatter(sensorLandTs, zscore(sensor(sensorLandTs, landSensor)), 20, 'r', 'filled')
                            plot(linspace(1,10e4, length(jerk)), zscore(jerk))
                            hold off
                            subplot(2,2,4)
                            
                            fprintf([num2str(size(jumpTs)) '\n'])
                            plot(p)
                            hold on
                            plot(jumpTs*ones(1,2), [min(p) max(p)], 'r') ;
                            hold off
                            shg
                            pause; %s{k}(i) = 'y';
                            prompt = input('Throw out: [Y]','s');
                            if isempty(prompt)
                                % convert to pos indices
                                maxImpLoc = jumpTs+jump_ts{k}(i)-120;
                                maxLandLoc = round((sensorLandTs-20e3)/167)+jump_ts{k}(i);
                                
                                xPos = movmean(pos(maxImpLoc-600:maxImpLoc+300,8), 5);
                                yPos = movmean(pos(maxImpLoc-600:maxImpLoc+300,10), 5);
                                b = 1;
                                
                                if sum(diff(xPos)) > 0
                                    mod = 1;
                                else
                                    mod = 2;
                                end
                                
                                
                                % Get start and end
                                if mod == 1 %
                                    startDist = pdist2([startx(1) starty(1)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                    endDist = pdist2([endx(2) endy(2)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                    [pks1, start_point] = findpeaks(-startDist, 'MinPeakHeight', -50, 'MinPeakWidth', 10);
                                    [pks2, end_point] = findpeaks(-endDist, 'MinPeakHeight', -50, 'MinPeakWidth', 10);
                                    
%                                     figure
%                                     subplot(3, 1, 1)
%                                     plot(-startDist)
%                                     hold on
%                                     scatter(start_point, pks1)
%                                     subplot(3, 1, 2)
%                                     plot(-endDist)
%                                     hold on
%                                     scatter(end_point, pks2)
%                                     subplot(3, 1, 3)
%                                     hold on
%                                     plot(xPos, yPos)
%                                     scatter(xPos(start_point), yPos(start_point))
%                                     scatter(xPos(end_point), yPos(end_point))
%                                     pause
%                                     close
                                    
                                    % Want it so that it chooses for
                                    % closest and for middle-est number
                                    %                                         [pks1, start_point] = min(abs(xPos(b:end-b)-endx(1))+abs(yPos(b:end-b)-endy(1)));
                                    %                                         [pks2, end_point] = min(abs(xPos(b:end-b)-endx(2))+ abs(yPos(b:end-b)-endy(2)));
                                    
                                    %                                         if pks1 > 50 || pks2 > 50
                                    %                                             [pks1, start_point] = min(pdist2([endx(1) endy(1)], [xPos yPos], 'euclidean'));
                                    %                                             [pks2, end_point] = min(pdist2([endx(2) endy(2)], [xPos yPos], 'euclidean'));
                                    %
                                    % %                                             [pks1, start_point] = min(abs(xPos-endx(1))+abs(yPos-endy(1)));
                                    % %                                             [pks2, end_point] = min(abs(xPos-endx(2))+ abs(yPos-endy(2)));
                                    %                                         end
                                else
                                    startDist = pdist2([startx(2) starty(2)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                    endDist = pdist2([endx(1) endy(1)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                    [pks1, start_point] = findpeaks(-startDist, 'MinPeakHeight', -50, 'MinPeakWidth', 10);
                                    [pks2, end_point] = findpeaks(-endDist, 'MinPeakHeight', -50, 'MinPeakWidth', 10);
                                    
                                    %                                         [pks2, end_point] = min(abs(xPos(b:end-b)-endx(1))+abs(yPos(b:end-b)-endy(1)));
                                    %                                         if pks1 > 50 || pks2 > 50
                                    %                                             [pks1, start_point] = min(pdist2([endx(2) endy(2)], [xPos yPos], 'euclidean'));
                                    %                                             [pks2, end_point] = min(pdist2([endx(1) endy(1)], [xPos yPos], 'euclidean'));
                                    %
                                    % %                                             [pks1, start_point] = min(abs(xPos-endx(2))+abs(yPos-endy(2)));
                                    % %                                             [pks2, end_point] = min(abs(xPos-endx(1))+ abs(yPos-endy(1)));
                                    %                                         end
                                end
                                
                                if ~isempty(start_point) && ~isempty(end_point)
                                    start_point = start_point(end)+b;
                                    end_point = end_point(1) +b;
                                    
                                    
                                    % Then check and store data
                                    % If z-range is below 500
                                    if max(pos(maxImpLoc-400:maxImpLoc+400,9)) - min(pos(maxImpLoc-400:maxImpLoc+400,9)) < 500
                                        behavior.events.trials{count}.x =pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,8);
                                        behavior.events.trials{count}.y =pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,10);
                                        behavior.events.trials{count}.z =pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,9);
                                        behavior.events.trials{count}.mapping = nan;
                                        behavior.events.trials{count}.timestamps = pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,1);
                                        behavior.events.trials{count}.errorPerMarker = pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,11);
                                        
                                        
                                        %             [r1 r2 r3] = quat2angle([pos(maxImpLoc-400:maxImpLoc+400,4);pos(maxImpLoc-400:maxImpLoc+400,5)...
                                        %                 pos(maxImpLoc-400:maxImpLoc+400,6);pos(maxImpLoc-400:maxImpLoc+400,7)])
                                        behavior.events.trials{count}.orientation.x = pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,4);
                                        behavior.events.trials{count}.orientation.y = pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,5);
                                        behavior.events.trials{count}.orientation.z = pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,6);
                                        behavior.events.trials{count}.orientation.w = pos(start_point+maxImpLoc-600:end_point+maxImpLoc-600,7);
                                        
                                        % for each jump, check which direction the rat ran
                                        %                                     if sum(diff(behavior.events.trials{count}.x)) > 0 %&& sum(diff(behavior.events.trials{count}.x)) > 0
                                        %                                         mod = 1;
                                        %                                     else
                                        %                                         mod = 2;
                                        %                                     end
                                        
                                        behavior.events.trialConditions(count) = k*2-(2-mod);% to make unique k %k * mod;
                                        behavior.events.trialIntervals(count,:) = pos([start_point+maxImpLoc-600 end_point+maxImpLoc-600],1);
                                        behavior.events.jumpTime(count,:) = [pos(maxImpLoc,1) pos(maxLandLoc,1)];
                                        behavior.events.jumpLoc(count, :) = [pos(maxImpLoc, 8) pos(maxLandLoc, 8) pos(maxImpLoc, 10) pos(maxLandLoc, 10)];
                                        count = 1+count;
                                    end
                                end
                            end
                            
                        end
                    end
                end
            end
        end
    end
    behavior.events.conditionType{k} = 'jump';
end
behavior.events.startEndPos = [startx starty endx endy];


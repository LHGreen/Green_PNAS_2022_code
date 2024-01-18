 function [behavior] = bz_getJumpBehav_piezo(pos, nchannels)

% parms = inputParser;
% addParameter(parms,'usePiezoSensor',false,@islogical);
% addParameter(parms,'nchannels',2,@isnumeric);
%
% parse(parms,varargin{:})

if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;
end
p = pos(:,[8 10]);
behavior = pos2behav(pos,'behavType','jump');

p(:,1) = movmean(p(:,1),40);
p(:,2) = movmean(p(:,2),40);

figure
scatter(pos(:,8),pos(:,9),'.')
% scatter(pos(:,8)-pos(:,10),pos(:,9),'.') % x-z is horizontal direction, y is height
axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
axis([-300 1200 500 900])%([-1000 2000 0 1000])
[x, y] = ginput(); % get click input


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
            if abs(pks_x{i}(j)) < 150 & abs(pks_y{i}(k)) < 20 % check that point is reasonably close
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
        if jump_ts{k}(i) > 1000 && jump_ts{k}(i) < size(pos, 1)-5000% !!!!!!!!!!!! Insert some likely condition here
            
            time = pos(jump_ts{k}(i),1);
            
            % get piezo sensor info
            sensor = bz_LoadBinary('analogin.dat', 'nchannels', nchannels, ...
                'start', time-1, 'duration', 5);
            
            sensor = (movmean(sensor, 40));
            %             [upper1, ~] = envelope(double(sensor(:,1)), 2200, 'peak');
            %             [upper2, ~] = envelope(double(sensor(:,2)), 2200, 'peak');
            
            % Find jump point
            %             [max1, id1] = max(upper1)
            %             [max2, id2] = max(upper2)
            [max1, id1] = max(abs(zscore(sensor(:,1))));
            [max2, id2] = max(abs(zscore(sensor(:,nchannels))));
            
            sensorJumpTs = 0;
            sensorMax = 0;
            
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
                    landTs = (round((sensorLandTs-20e3)/167)+jump_ts{k}(i));
                    
                    jump_ts{k}(i)
                    p = pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,9); % position
                    px = pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,8); % position
                    jerk = diff(diff(fastrms(fastrms(diff(movmean(px,40)), 40),80)));
                    
                    [~, jumpTs] = max(jerk);
                    
                    %             if abs(id1-id2) < 10000 && abs(id1-id2) > 2000 % must be within 100 and 500 ms
                    %                 if id1 < id2
                    %                     if max1 > mean(sensor(:,1))+2*sqrt(var(sensor(:,1)))
                    %                         sensorJumpTs = id1;
                    %                         sensorMax = max1;
                    %                         sensorLandTs = id2;
                    %                         sensorMax2 = max2;
                    %                     end
                    %                 elseif id1 > id2
                    %                     if max2 > mean(sensor(:,1))+2*sqrt(var(sensor(:,1)))
                    %                         sensorJumpTs = id2;
                    %                         sensorMax = max2;
                    %                         sensorLandTs = id1;
                    %                         sensorMax2 = max1;
                    %                     end
                    %                 end
                    %             end
                    
                    %             if sensorMax ~=0
                    if (jumpTs+jump_ts{k}(i)-120) < landTs
                        if abs((jumpTs+jump_ts{k}(i)-120)-landTs) < 200 && abs((jumpTs+jump_ts{k}(i)-120)-landTs) > 30
                            
                            
                            % Start by plotting the putative points
                            % plot position in horizontal plane
                            subplot(2,2,1)
                            scatter(pos(:,8),pos(:,10),'.')
                            hold on
                            % plot position around jump time only
                            scatter(pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,8),pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,10),'.r')
                            title(['count = ', num2str(count)])
                            hold off
                            subplot(2,2,2)
                            % plot position in the xy direction
                            scatter(pos(:,8),pos(:,9),'.')
                            hold on
                            % plot position around the jump time only
                            scatter(pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,8),pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,9),'.r')
                            axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
                            hold off
                            % plot sensor input
                            %                     p = pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,9); % position
                            %                     px = pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,8); % position
                            %                     jerk = diff(diff(fastrms(fastrms(diff(movmean(px,40)), 40),80)));
                            subplot(2,2,3)
                            plot(zscore(sensor(:,1)))
                            hold on
                            plot(zscore(sensor(:,2)))
                            %                     plot(upper1, 'k', 'LineWidth', 0.25)
                            %                     plot(upper2, 'k', 'LineWidth', 0.25)
                            scatter(jumpTs*167, zscore(jerk(jumpTs)), 20, 'g', 'filled')
                            scatter(sensorLandTs, zscore(sensor(sensorLandTs, landSensor)), 20, 'r', 'filled')
                            plot(linspace(1,10e4, length(jerk)), zscore(jerk))
                            hold off
                            subplot(2,2,4)
                            
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
                                % Then check and store data
                                if max(pos(maxImpLoc-400:maxImpLoc+400,9)) - min(pos(maxImpLoc-400:maxImpLoc+400,9)) < 500
                                    behavior.events.trials{count}.x =pos(maxImpLoc-400:maxImpLoc+400,8);
                                    behavior.events.trials{count}.y =pos(maxImpLoc-400:maxImpLoc+400,10);
                                    behavior.events.trials{count}.z =pos(maxImpLoc-400:maxImpLoc+400,9);
                                    behavior.events.trials{count}.mapping = 1:801;
                                    behavior.events.trials{count}.timestamps = pos(maxImpLoc-400:maxImpLoc+400,1);
                                    behavior.events.trials{count}.errorPerMarker = pos(maxImpLoc-400:maxImpLoc+400,11);
                                    
                                    
                                    %             [r1 r2 r3] = quat2angle([pos(maxImpLoc-400:maxImpLoc+400,4);pos(maxImpLoc-400:maxImpLoc+400,5)...
                                    %                 pos(maxImpLoc-400:maxImpLoc+400,6);pos(maxImpLoc-400:maxImpLoc+400,7)])
                                    behavior.events.trials{count}.orientation.x = pos(maxImpLoc-400:maxImpLoc+400,4);
                                    behavior.events.trials{count}.orientation.y = pos(maxImpLoc-400:maxImpLoc+400,5);
                                    behavior.events.trials{count}.orientation.z = pos(maxImpLoc-400:maxImpLoc+400,6);
                                    behavior.events.trials{count}.orientation.w = pos(maxImpLoc-400:maxImpLoc+400,7);
                                    
                                    % for each jump, check which direction the rat ran
                                    if sum(diff(behavior.events.trials{count}.x)) > 0 & sum(diff(behavior.events.trials{count}.x)) > 0
                                        mod = 1;
                                    else
                                        mod = 2;
                                    end
                                    
                                    behavior.events.trialConditions(count) = k*2-(2-mod);% to make unique k %k * mod;
                                    behavior.events.trialIntervals(count,:) = pos([maxImpLoc-400 maxImpLoc+400],1);
                                    behavior.events.jumpTime(count,:) = [pos(maxImpLoc,1) pos(maxLandLoc,1)];
                                    count = 1+count;
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














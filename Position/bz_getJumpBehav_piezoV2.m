function [behavior] = bz_getJumpBehav_piezoV2(pos, nchannels)

% I'm not sure what the plots show anymore.
% parms = inputParser;
% addParameter(parms,'usePiezoSensor',false,@islogical);
% addParameter(parms,'nchannels',2,@isnumeric);
%
% parse(parms,varargin{:})

if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;
end
p = pos(:,[8 9]); % x y
behavior = pos2behav(pos,'behavType','jump');

p(:,1) = movmean(p(:,1),10);
p(:,2) = movmean(p(:,2),10);

% %% Use Kalman filter to smooth the x position
% % Kalman filter
% % Constant Velocity Model
% % Using Adaptive Adjustment of Noise Covariance in Kalman Filter for
% % Dynamic Estimation - Akhlagi, Zhou, & Huang, Arxiv
% 
% % g = 9.18*10e3/120; % 9.18 m/s
% T = (nanmean(diff(pos(:,1)))); % timestep, should be 1/120
% 
% u = 0; %-g;
% x0 = [pos(1,8) diff(pos(1:2, 8))]'; % x-position start
% A = [1 T ; 0 1]; %
% B = [ T 1]';
% 
% C = [1 0];
% 
% error = pos(:, 11);
% error(isnan(error)) = [];
% 
% R = cov(error); % initial estimate
% alpha = 0.3;
% 
% C = [1 0]; % output transition matrix
% 
% Q = 1*eye(2); % process noise covariance
% R = 53; % measurement noise covariance
% 
% 
% 
% % Implement the filter
% t = length(pos(:,1));
% ye = zeros(t, 1);
% ycov = zeros(t, 1);
% yv = pos(1:t,8); % the z-position
% x = x0;
% P = cov([yv(1:2) zeros(2, 1)]); % covariance matrix
% 
% %
% for ii = 1:t
%     % measurement update
%     if sum(isnan(x)) > 0
%         x = zeros(2, 1);
%     end
%     
%     % Step 1: Prediction update
%     x = A*x + B*u;
%     P = A*P*A' + Q;
%     
%     % Step 2: Correction
%     d = yv(ii)-C*x; % measurement innovation
%     S = C*P*C' + R; % Innovation Covariance
%     K = P*C'/S; % Kalman gain matrix
%     x = x + K*d;
%     P_m = P; % save for next step
%     P = (eye(2) -K*C)*P; %*(eye(2)-Mn*C)' + Mn*R*Mn';
%     
%     ye(ii) = C*x;
%     
%     ycov(ii) = C*P*C';
%     
%     % estimate R and Q
%     resid = yv(ii)-C*x;
%     S_hat = resid*resid';
%     %     R = alpha*R + (1-alpha)*(S_hat + C*P_m*C');
%     
%     Q = alpha*Q + (1-alpha)*(K*(d*d')*K');
%     
%     % EM algorithm to estimate Q
%     
%     %     X = P + x*x';
%     %     J = P*A'*(A*P*A'+Q);
%     %     XT =  J*(xx-A*x);
%     
%     
%     
%     
% end
% %
% figure
% plot(yv)
% hold on
% plot(ye)
% pause
% close
% 
% p(:,1) = ye; % replace position vector
% 


%%


figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(p(:,1),p(:,2),'.') % x-z is horizontal direction, y is height
% axis([min(pos(:,8)) max(pos(:,8)) min(pos(:,9)) max(pos(:,9))])
axis([-500 600 400 1200])%([-1000 2000 0 1000])
[x, y] = ginput(); % get click input


% Then get start/stop positions
figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(p(:,1),pos(:,10),'.') % x-z is horizontal direction, y is height
axis([-800 1400 -1400 1200])%
% axis([600 2800 -1600 1200])%
title('start points (L, then R)')
[startx, starty] = ginput(2);
close

figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(p(:,1),pos(:,10),'.') % x-z is horizontal direction, y is height
axis([-800 1400 -1400 1200])%
% axis([600 2800 -1600 1200])
title('end points (L, then R)')
[endx, endy] = ginput(2);

close

for i=1:length(x) % for all mouse clicks
    dists_x(i,:) = abs(p(:,1)-x(i));
    dists_y(i,:) = abs(pos(:,9)-y(i));
    
    % find points in pos that are closest to mouse click
    [pks_x{i} locs_x{i}] = findpeaks(-dists_x(i,:),'MinPeakHeight',-300);
    [pks_y{i} locs_y{i}] = findpeaks(-dists_y(i,:),'MinPeakHeight',-70);
    for j=1:length(x)
        trials{i,j,1}=[];
        trials{i,j,2}=[];
    end
    jump_ts{i} = 0;
    % for every peak that is close to a mouse click
    for j=length(pks_x{i}):-1:1
        for k=length(pks_y{i}):-1:1
            if abs(pks_x{i}(j)) < 60 && abs(pks_y{i}(k)) < 60 % check that point is reasonably close
                if abs(locs_x{i}(j) - locs_y{i}(k)) < 20 % if aligned in time
                    if abs(mean([locs_x{i}(j) locs_y{i}(k)]) - jump_ts{i}(end)) > 180 % if at least 1.5 s after previous jump_ts
                        if size(pos, 1) - mean([locs_x{i}(j) locs_y{i}(k)]) > 480 % if not at end of trial
                            jump_ts{i} = [jump_ts{i};round(mean([locs_x{i}(j) locs_y{i}(k)]))]; % jump ts is mean of x and y indices
                        end
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
        if jump_ts{k}(i) > 121 % !!!!!!!!!!!! Insert some likely condition here
            
            time = pos(jump_ts{k}(i),1);
            
            % get piezo sensor info
            sensor = bz_LoadBinary('analogin.dat', 'nchannels', nchannels, ...
                'start', time-1, 'duration', 3.5);
            
            sensor = (movmean(sensor, 40));
            %
            [max1, id1] = max(abs(zscore(sensor(:,1))));
            [max2, id2] = max(abs(zscore(sensor(:,2))));
            
            sensorJumpTs = 0;
            sensorMax = 0;
            
            % Determine which sensor recorded the landing
            landSensor = 0;
            if id1 < id2
                landSensor = 2;
            elseif id1 > id2
                landSensor = 1;
            end
            
            if landSensor ~=0
                if jump_ts{k}(i) < size(pos, 1) -1000 && jump_ts{k}(i) > 500
                    meanSensor = mean(sensor(8e3:10e3,landSensor));
                    sdSensor = sqrt(var(sensor(8e3:10e3,landSensor)));
                    sensorLandTs = find(abs(sensor(:,landSensor)-meanSensor) > 7*sdSensor,1);
                    landTs = (round((sensorLandTs-20e3)/167)+jump_ts{k}(i));
                    
                    jump_ts{k}(i)
                    pz = pos(jump_ts{k}(i)-120:jump_ts{k}(i)+480,9); % z-position
                    px = p(jump_ts{k}(i)-120:jump_ts{k}(i)+480,1); % horizontal position
                    vel = smoothdata(diff(smoothdata(px, 'sgolay', 5)),'sgolay', 5);
                    acc = smoothdata(diff(vel), 'sgolay', 5);
                    jerk = diff(diff(fastrms(fastrms(diff(movmean(px,40)), 40),80)));
                    
                    %                     [~, jumpTs] = max(jerk);
                    %                     [~, jumpTs] = max(abs(acc));
                    
                    t = round(sensorLandTs/167);
                    
                    jumpTs = [];
                    
                    if t > 60
                        
                        [~, jumpTs] = max(abs(acc((t-60):t)));
                        jumpTs = jumpTs + t-61;
                    elseif t > 0
                        [~, jumpTs] = max(abs(acc(1:t)));
                        jumpTs = jumpTs + t-1;
                    end
                    
                    if isempty(jumpTs) || isempty(landTs)
                        close
                        
                    else
                        
                        
                        if (jumpTs+jump_ts{k}(i)-120) < landTs
                            % check that jump duration is less than 1 s and longer than 0.100 seconds 
                            if abs((jumpTs+jump_ts{k}(i)-120)-landTs) < 1*120  && abs((jumpTs+jump_ts{k}(i)-120)-landTs) > 0.10*120
                                
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
                                axis([min(pos(:,8)) max(pos(:,8)) min(pos(:,9)) max(pos(:,9))])
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
                                scatter(jumpTs*167, zscore(acc(jumpTs)), 40, 'g', 'filled')
                                scatter(sensorLandTs, zscore(sensor(sensorLandTs, landSensor)), 40, 'r', 'filled')
                                plot(linspace(1,10e4, length(acc)), zscore(acc))
                                hold off
                                legend('sensor 1', 'sensor 2', 'jump', 'land', 'acc')
                                subplot(2,2,4)
                                plot(px)
                                hold on
                                plot(jumpTs*ones(1,2), [min(px) max(px)], 'r') ;
                                plot(sensorLandTs*ones(1, 2)*(1/167), [min(px) max(px)], 'g') ;
                                hold off
                                shg
                                pause; %s{k}(i) = 'y';
                                
                                
                                prompt = input('Throw out: [Y]','s');
                                if isempty(prompt)
                                    % convert to pos indices
                                    maxImpLoc = jumpTs+jump_ts{k}(i)-120;
                                    maxLandLoc = round((sensorLandTs-20e3)/167)+jump_ts{k}(i);
                                    
                                    % Then refine trial positions
                                    buffer = 600;
                                    % use these values when the rat pauses a long
                                    % time
                                    xPos = movmean(pos(maxImpLoc-buffer:maxImpLoc+300,8), 5);
                                    yPos = movmean(pos(maxImpLoc-buffer:maxImpLoc+300,10), 5);
                                    b = 1;
                                    
                                    if sum(diff(xPos)) > 0
                                        mod = 1;
                                    else
                                        mod = 2;
                                    end
                                    
                                    % Get start and end
                                    if mod == 1 %
                                        % get start points closest to the start
                                        % points
                                        startDist = pdist2([startx(1) starty(1)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                        endDist = pdist2([endx(2) endy(2)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                        [pks1, start_point] = findpeaks(-startDist, 'MinPeakHeight', -60, 'MinPeakWidth', 10);
                                        [pks2, end_point] = findpeaks(-endDist, 'MinPeakHeight', -60, 'MinPeakWidth', 10);
                                        
                                        
                                    else
                                        startDist = pdist2([startx(2) starty(2)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                        endDist = pdist2([endx(1) endy(1)], [xPos(b:end-b) yPos(b:end-b)], 'euclidean');
                                        [pks1, start_point] = findpeaks(-startDist, 'MinPeakHeight', -60, 'MinPeakWidth', 10);
                                        [pks2, end_point] = findpeaks(-endDist, 'MinPeakHeight', -60, 'MinPeakWidth', 10);
                                        
                                        
                                    end
                                    
                                    if ~isempty(start_point) && ~isempty(end_point)
                                        start_point = start_point(end)+b;
                                        end_point = end_point(1) +b;
                                        
                                        % Check that points are within bounds
                                        if min(pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer, 8)) > (min(endx)-30) & max(pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer, 8)) < (max(endx)+30)
                                            
                                            % Then check and store data
                                            % If z-range is below 500
                                            if max(pos(maxImpLoc-400:maxImpLoc+400,9)) - min(pos(maxImpLoc-400:maxImpLoc+400,9)) < 500
                                                behavior.events.trials{count}.x =pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,8);
                                                behavior.events.trials{count}.y =pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,10);
                                                behavior.events.trials{count}.z =pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,9);
                                                behavior.events.trials{count}.mapping = nan;
                                                behavior.events.trials{count}.timestamps = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,1);
                                                behavior.events.trials{count}.errorPerMarker = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,11);
                                                
                                                behavior.events.trials{count}.orientation.x = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,4);
                                                behavior.events.trials{count}.orientation.y = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,5);
                                                behavior.events.trials{count}.orientation.z = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,6);
                                                behavior.events.trials{count}.orientation.w = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,7);
                                                
                                                behavior.events.trialConditions(count) = k*2-(2-mod);% to make unique k %k * mod;
                                                behavior.events.trialIntervals(count,:) = pos([start_point+maxImpLoc-buffer end_point+maxImpLoc-buffer],1);
                                                behavior.events.jumpTime(count,:) = [pos(maxImpLoc,1) pos(maxLandLoc,1)];
                                                behavior.events.jumpLoc(count, :) = [pos(maxImpLoc, 8) pos(maxLandLoc, 8) pos(maxImpLoc, 10) pos(maxLandLoc, 10) pos(maxImpLoc, 9) pos(maxLandLoc, 9)];
                                                count = 1+count;
                                            end
                                        end
                                    end
                                    
                                    
                                    
                                    
                                    
                                    %                                 % Then check and store data
                                    %                                 if max(pos(maxImpLoc-400:maxImpLoc+400,9)) - min(pos(maxImpLoc-400:maxImpLoc+400,9)) < 500
                                    %                                     behavior.events.trials{count}.x =pos(maxImpLoc-400:maxImpLoc+400,8);
                                    %                                     behavior.events.trials{count}.y =pos(maxImpLoc-400:maxImpLoc+400,10);
                                    %                                     behavior.events.trials{count}.z =pos(maxImpLoc-400:maxImpLoc+400,9);
                                    %                                     behavior.events.trials{count}.mapping = 1:801;
                                    %                                     behavior.events.trials{count}.timestamps = pos(maxImpLoc-400:maxImpLoc+400,1);
                                    %                                     behavior.events.trials{count}.errorPerMarker = pos(maxImpLoc-400:maxImpLoc+400,11);
                                    %
                                    %
                                    %                                     %             [r1 r2 r3] = quat2angle([pos(maxImpLoc-400:maxImpLoc+400,4);pos(maxImpLoc-400:maxImpLoc+400,5)...
                                    %                                     %                 pos(maxImpLoc-400:maxImpLoc+400,6);pos(maxImpLoc-400:maxImpLoc+400,7)])
                                    %                                     behavior.events.trials{count}.orientation.x = pos(maxImpLoc-400:maxImpLoc+400,4);
                                    %                                     behavior.events.trials{count}.orientation.y = pos(maxImpLoc-400:maxImpLoc+400,5);
                                    %                                     behavior.events.trials{count}.orientation.z = pos(maxImpLoc-400:maxImpLoc+400,6);
                                    %                                     behavior.events.trials{count}.orientation.w = pos(maxImpLoc-400:maxImpLoc+400,7);
                                    %
                                    %                                     % for each jump, check which direction the rat ran
                                    %                                     if sum(diff(behavior.events.trials{count}.x)) > 0 & sum(diff(behavior.events.trials{count}.x)) > 0
                                    %                                         mod = 1;
                                    %                                     else
                                    %                                         mod = 2;
                                    %                                     end
                                    %
                                    %                                     behavior.events.trialConditions(count) = k*2-(2-mod);% to make unique k %k * mod;
                                    %                                     behavior.events.trialIntervals(count,:) = pos([maxImpLoc-400 maxImpLoc+400],1);
                                    %                                     behavior.events.jumpTime(count,:) = [pos(maxImpLoc,1) pos(maxLandLoc,1)];
                                    %                                     count = 1+count;
                                    %                                 end
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



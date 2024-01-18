function [behavior] = bz_getJumpBehav_NoPiezoEditVersion(pos)

dbstop if error

if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;
end
behavior = pos2behav(pos,'behavType','jump');

% pos(:,8:10) = pos(:,8:10)*1000;

%% Smooth z-position
% % Kalman filter
% % Constant Velocity Model
% % Using Adaptive Adjustment of Noise Covariance in Kalman Filter for
% % Dynamic Estimation - Akhlagi, Zhou, & Huang, Arxiv
% 
% % g = 9.18*10e3/120; % 9.18 m/s
T = (nanmean(diff(pos(:,1)))); % timestep, should be 1/120
% 
% u = 0; %-g;
% 
% x0 = [pos(1,9) 0]'; % z-position start
% A = [1 T ; 0 1]; % Transition Matrix
% B = [ T 1]'; % input matrix (maps control commands onto state changes)
% % % equation x = Ax + Bu + w, u is input control vector, w is noise
% % % assume observations are z = Cx + v, where C is observation matrix, v is
% % % measurement noise
% % % because w has mean 0, the predictor is x = Ax + Bu.
% % % covariance matrix P = APA' + Q
% % 
% % error = pos(:, 11);
% % error(isnan(error)) = [];
% % R = cov(error); % initial estimate
% alpha = 0.99;

% C = [1 0]; % output transition matrix
% 
% Q = 1*eye(2); % process noise covariance
% R = 15; %30; % measurement noise covariance
% 
% % Implement the filter
% t = length(pos(:,1));
% 
% 
% ye = zeros(t, 1);
% ycov = zeros(t, 1);
% yv = pos(1:t,9); % the z-position
% x = x0; 
% P = eye(2);% covariance matrix
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
% %     R = alpha*R + (1-alpha)*(S_hat + C*P_m*C');
%     
%     Q = alpha*Q + (1-alpha)*(K*(d*d')*K');
%     
%     % EM algorithm to estimate Q
%     
% %     X = P + x*x';
% %     J = P*A'*(A*P*A'+Q);
% %     XT =  J*(xx-A*x);
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
% posZ = ye; % replace position
% 
% %% Do the same for x-position
% 
% x0 = [pos(1,8) diff(pos(1:2, 8))]'; % x-position start
% A = [1 T ; 0 1]; % 
% B = [ T 1]';
% 
% 
% % R = cov(error); % initial estimate
% alpha = 0.8;
% 
% C = [1 0]; % output transition matrix
% 
% Q = 1*eye(2); % process noise covariance
% R = 15; % measurement noise covariance
% 
% 
% 
% % Implement the filter
% t = length(pos(:,1));
% ye = zeros(t, 1);
% ycov = zeros(t, 1);
% yv = pos(1:t,8); % the x-position
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
% %     R = alpha*R + (1-alpha)*(S_hat + C*P_m*C');
%     
%     Q = alpha*Q + (1-alpha)*(K*(d*d')*K');
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
% posX = ye; % replace position
% 
% %% Same for y pos now
% 
% 
% x0 = [pos(1,10) diff(pos(1:2, 10))]'; % x-position start
% A = [1 T ; 0 1]; % 
% B = [ T 1]';
% 
% 
% % R = cov(error); % initial estimate
% alpha = 0.8;
% 
% C = [1 0]; % output transition matrix
% 
% Q = 1*eye(2); % process noise covariance
% R = 15; % measurement noise covariance
% 
% 
% 
% % Implement the filter
% t = length(pos(:,1));
% ye = zeros(t, 1);
% ycov = zeros(t, 1);
% yv = pos(1:t,10); % the y-position
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
% %     R = alpha*R + (1-alpha)*(S_hat + C*P_m*C');
%     
%     Q = alpha*Q + (1-alpha)*(K*(d*d')*K');
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
% posY = ye; % replace position


%% Smooth position
% Kalman filter
% Constant Velocity Model
% Using Adaptive Adjustment of Noise Covariance in Kalman Filter for
% Dynamic Estimation - Akhlagi, Zhou, & Huang, Arxiv

% It's probably better to smooth all of them?
% 
% T = (nanmean(diff(pos(:,1)))); % timestep, should be 1/120
% 
% u = 0; %-g;
% 
% x0 = [pos(1,8) 0 pos(1, 10), 0 pos(1, 9), 0]'; % X, Y, Z-position start
% A = [1 T 0 0 0 0; 0 1 0 0 0 0; ...
%     0 0 1 T 0 0; 0 0 0 1 0 0; ...
%     0 0 0 0 1 T; 0 0 0 0 0 1]; % Transition Matrix
% B = [ T 1 0 0 0 0]'; % input matrix (maps control commands onto state changes)
% % equation x = Ax + Bu + w, u is input control vector, w is noise
% % assume observations are z = Cx + v, where C is observation matrix, v is
% % measurement noise
% % because w has mean 0, the predictor is x = Ax + Bu.
% % covariance matrix P = APA' + Q
% 
% error = pos(:, 11);
% error(isnan(error)) = [];
% R = cov(error); % initial estimate
% alpha = 0.3;
% 
% C = [1 0 0 0 0 0; ...
%     0 0 1 0 0 0; ...
%     0 0 0 0 1 0];
% 
% Q = 1*eye(6); % process noise covariance
% R = 47; % measurement noise covariance
% 
% % Implement the filter
% t = length(pos(:,1));
% 
% 
% ye = zeros(t, 3);
% ycov = zeros(t, 1);
% vv = [0 0 0; diff(pos(1:t, 8:10))];
% yv = [pos(1:t,8) vv(:, 1) pos(1:t, 10), vv(:, 3), pos(1:t, 9) vv(:, 2)]; % the z-position
% x = x0;
% P = cov(yv(1:2, :)); % covariance matrix
% 
% %
% for ii = 1:t
%     % measurement update
%     if sum(isnan(x)) > 0
%         x = zeros(6, 1);
%     end
%     
%     % Step 1: Prediction update
%     x = A*x + B*u;
%     P = A*P*A' + Q;
%     
%     % Step 2: Correction
%     d = yv(ii,[1 3 5])'-C*x; % measurement innovation
%     S = C*P*C' + R; % Innovation Covariance
%     K = P*C'*inv(S); % Kalman gain matrix
%     x = x + K*d;
%     P_m = P; % save for next step
%     P = (eye(6) -K*C)*P; %*(eye(2)-Mn*C)' + Mn*R*Mn';
%     
%     out = C*x;
%     ye(ii, :) = out;
%     
% %     ycov(ii) = C*P*C';
%     
%     % estimate R and Q
%     resid = yv(ii, [1 3 5])-C*x;
%     S_hat = resid*resid';
% %     R = alpha*R + (1-alpha)*(S_hat + C*P_m*C');
%     
%     Q = alpha*Q + (1-alpha)*(K*(d*d')*K');
%     
%     % EM algorithm to estimate Q
%     
% %     X = P + x*x';
% %     J = P*A'*(A*P*A'+Q);
% %     XT =  J*(xx-A*x);
%     
%     
%     
%     
% end
% %
% figure
% plot(yv(:, 1))
% hold on
% plot(ye(:, 1))
% % pause
% % close
% 
% figure
% plot(yv(:, 2))
% hold on
% plot(ye(:, 2))
% 
% figure
% plot(yv(:, 3))
% hold on
% plot(ye(:, 3))
% pos(:,9) = ye; % replace position


%% Kalman smooth Z position

x0 = [pos(1,9) diff(pos(1:2, 9))]'; % x-position start
A = [1 T ; 0 1]; % 
B = [ T 1]';

u = 0;
C = [1 0];


Q = 90; % process noise covariance
R = 1; % measurement noise covariance

% Implement the filter
t = length(pos(:,1));

P = B*Q*B'; % covariance matrix
ye = zeros(t, 1);
ycov = zeros(t, 1);
yv = pos(1:t,9); % the z-position
x = x0;


for ii = 1:t
    % measurement update
    if sum(isnan(x)) > 0
        x = zeros(2, 1);
    end
    Mn = P*C'/(C*P*C'+R); % Kalman gain matrix
    x = x + Mn*(yv(ii)-C*x);
    P = (eye(2) -Mn*C)*P;
    
    ye(ii) = C*x;
    
    ycov(ii) = C*P*C';
    
    % Time update
    x = A*x + B*u;
    P = A*P*A' + B*Q*B';
end
%
figure
plot(yv)
hold on
plot(ye)
pause
close

posZ = ye; % replace position

%% Do the same for x-position

x0 = [pos(1,8) diff(pos(1:2, 8))]'; % x-position start
A = [1 T ; 0 1]; % 
B = [ T 1]';

C = [1 0];

Q = 25;
R = 20;

% Implement the filter
t = length(pos(:,1));

P = B*Q*B';
ye = zeros(t, 1);
ycov = zeros(t, 1);
yv = pos(1:t,8); % the x-position
x = x0;


for ii = 1:t
    % measurement update
    if sum(isnan(x)) > 0
        x = zeros(2, 1);
    end
    Mn = P*C'/(C*P*C'+R);
    x = x + Mn*(yv(ii)-C*x);
    P = (eye(2) -Mn*C)*P;
    
    ye(ii) = C*x;
    
    ycov(ii) = C*P*C';
    
    % Time update
    x = A*x + B*u;
    P = A*P*A' + B*Q*B';
end
%
figure
plot(yv)
hold on
plot(ye)
pause
close

posX = ye; % replace position


%% Do the same for y-position

x0 = [pos(1,10) diff(pos(1:2, 10))]'; % y-position start
A = [1 T ; 0 1]; % 
B = [ T 1]';

C = [1 0];

Q = 25;
R = 20;

% Implement the filter
t = length(pos(:,1));

P = B*Q*B';
ye = zeros(t, 1);
ycov = zeros(t, 1);
yv = pos(1:t,10); % the y-position
x = x0;


for ii = 1:t
    % measurement update
    if sum(isnan(x)) > 0
        x = zeros(2, 1);
    end
    Mn = P*C'/(C*P*C'+R);
    x = x + Mn*(yv(ii)-C*x);
    P = (eye(2) -Mn*C)*P;
    
    ye(ii) = C*x;
    
    ycov(ii) = C*P*C';
    
    % Time update
    x = A*x + B*u;
    P = A*P*A' + B*Q*B';
end
%
figure
plot(yv)
hold on
plot(ye)
pause
close

posY = ye; % replace position





%% 
% Then get start/stop positions
figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(posX,posY,'.') % x-z is horizontal direction, y is height
axis([-800 1400 -1400 1200])%
% axis([600 2800 -1600 1200])%
title('start points')
[startx, starty] = ginput(2);
close


figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(posX,posY,'.') % x-z is horizontal direction, y is height
axis([-800 1400 -1400 1200])%
% axis([600 2800 -1600 1200])
title('end points')
[endx, endy] = ginput(2);

close

%% get horizontal position

line = [diff(startx) diff(starty)]';

position = [posX posY]*(line/norm(line));

%%

% get jump positions
figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(position,posZ, '.')
% scatter3(pos(:,8), pos(:,9), pos(:,10), '.')
% scatter3(pos(:,8)+pos(:,10),pos(:,9),'.') % x-z is horizontal direction, y is height
% axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
axis([-800 1000 400 1000])%
% axis([400 1000 500 800])
% axis([2000 3200 400 1000])
% axis([800 1800 400 1000])
[x, y] = ginput(); % get click input
close



for i=1:length(x) % for all mouse clicks
    %     dists_x(i,:) = abs(pos(:,8)+pos(:,10)-x(i));
    %     dists_y(i,:) = abs(pos(:,8)+pos(:,10)-y(i));
    dists_x(i,:) = abs(position-x(i));
    dists_y(i,:) = abs(posZ-y(i));
    
    % find points in pos that are closest to mouse click
    [pks_x{i} locs_x{i}] = findpeaks(-dists_x(i,:),'MinPeakHeight',-300);
    [pks_y{i} locs_y{i}] = findpeaks(-dists_y(i,:),'MinPeakHeight',-70);
    %     for j=1:length(x)
    %         trials{i,j,1}=[];
    %         trials{i,j,2}=[];
    %     end
    jump_ts{i} = 0;
    % for every peak that is close to a mouse click
    for j=length(pks_x{i}):-1:1
        for k=length(pks_y{i}):-1:1
            if abs(pks_x{i}(j)) < 60 && abs(pks_y{i}(k)) < 60 % check that point is reasonably close to mouse click
                if abs(locs_x{i}(j) - locs_y{i}(k)) < 20 % if aligned in time
                    if abs(mean([locs_x{i}(j) locs_y{i}(k)]) - jump_ts{i}(end)) > 240 % if at least 240 units (~2 s) after previous jump_ts
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
        if jump_ts{k}(i) > 121
            
            % Instead of using sensor, just pick a spot
            
            pz = round(posZ(jump_ts{k}(i)-120:jump_ts{k}(i)+480), 5); % z-position
            px = round(position(jump_ts{k}(i)-120:jump_ts{k}(i)+480),5); % x-position
            figure('units', 'normalized', 'outerposition', [0 0 1 1])
            plot(pz)
            [t, z] = ginput(1);
            t = round(t);
            
            
            
            
            if ~isempty(t)
                if jump_ts{k}(i) < size(pos, 1) -1000
                    % this is landTs in terms of behavioral timestamps
                
                    vel = smoothdata(diff(smoothdata(px, 'gaussian', [0 3])),'gaussian', [0 3]);
                    acc = smoothdata(diff(vel), 'gaussian', [0 3]);
                    
                    
                    %                     rot = r(jump_ts{k}(i)-120:jump_ts{k}(i)+480);
                    %                     rval = r2(jump_ts{k}(i)-120:jump_ts{k}(i)+480);
                    
                    % z-position variables
                    pz = smoothdata(pz, 'gaussian', [0 5]); 
                    velz = smoothdata(diff(pz),'gaussian', [0 5]);
                    accz = smoothdata(diff(velz), 'gaussian', [0 5]);
                    acc_zscore = zscore(accz);
                    
                   
                    fprintf(['t = ' num2str(t)])
                    
                    if t > 60
                        
                        [~, jumpTs] = max(abs(acc((t-60):t)));
                        [~, landTs] = findpeaks(acc_zscore(t:-1:t-30), 'NPeaks', 1, 'MinPeakHeight', 0.4);
                        landTs = t-landTs+1;
                        jumpTs = jumpTs + t-61;
                    elseif t > 0
                        [~, jumpTs] = max(abs(acc(1:t)));
                        [~, landTs] = findpeaks(acc_zscore(t:-1:1), 'NPeaks', 1, 'MinPeakHeight', 0.4);
                        landTs = t-landTs+1;
                        jumpTs = jumpTs + t-1;
                    end
                    
                    if isempty(jumpTs) || isempty(landTs)
                        close
                        
                    else
                        %                     if (jumpTs+jump_ts{k}(i)-120) < landTs % if jump time is before land time
                        % if jump time is at least 0.1 ms and no greater
                        % than 0.5 ms
                        %                         if abs((jumpTs+jump_ts{k}(i)-120)-landTs) < 200 && abs((jumpTs+jump_ts{k}(i)-120)-landTs) > 12
                        
                        % Start by plotting the putative points
                        % plot position in horizontal plane
                        
                        
                        %% This is just a graphics fig
                        inds = jumpTs-50:jumpTs+50;
                        subplot(2,2,2)
                        plot(1:101, zscore(pz(inds)))
                        hold on
                        plot(1:101, zscore(velz(inds)))
                        plot(1:101, zscore(accz(inds)))
                        scatter(51, zscore(acc(jumpTs)), 20, 'g', 'filled')
                        %back to normal
                        scatter(landTs-jumpTs + 51, zscore(accz(landTs)), 20, 'r', 'filled')
                        box off
                        xlim([0 121])
                        legend('pz', 'velz', 'accz', 'jump')
                        title('z')
                        hold off
                        sgtitle(['jump # = ' num2str(k) ', trial = ' num2str(i) '/' num2str(length(jump_ts{k}))])
                        
                        subplot(2,2,1)
                       plot(1:101, zscore(px(inds)))
                       hold on
                        plot(1:101, zscore(vel(inds)))
                        plot(1:101, zscore(acc(inds)))
                        scatter(51, zscore(acc(jumpTs)), 20, 'g', 'filled')
                        scatter(landTs-jumpTs+ 51, zscore(acc(landTs)), 20, 'r', 'filled')
                        box off
                        xlim([0 121])
                        %                             plot(linspace(1, 10e4, length(jerk2)), jerk2/max(jerk2))
                        legend('jumpTime', 'landTime', 'px', 'vel', 'acc')
                        scatter(t, zscore(accz(t)), 20, 'b', 'filled')
                        hold off
                        title('x')
                        
                        pause
                        
                        close
                        
                        %%
                        
                        
                        subplot(2,2,2)
                        plot(1:101, zscore(pz(landTs-50:landTs+50)))
                        hold on
                        plot(1:101, zscore(velz(landTs-50:landTs+50)))
                        plot(1:101, zscore(accz(landTs-50:landTs+50)))
                        scatter((jumpTs-landTs)+51, zscore(acc(jumpTs)), 20, 'g', 'filled')
                        %back to normal
                        scatter(51, zscore(accz(landTs)), 20, 'r', 'filled')
                        scatter((t-landTs)+51, zscore(accz(t)), 20, 'b', 'filled')
                        legend('pz', 'velz', 'accz', 'jump', 'land')
                        title('z')
                        hold off
                        sgtitle(['jump # = ' num2str(k) ', trial = ' num2str(i) '/' num2str(length(jump_ts{k}))])
                        
                        subplot(2,2,1)
                       
                        scatter(jumpTs, zscore(acc(jumpTs)), 20, 'g', 'filled')
                        hold on
                        scatter(landTs, zscore(acc(landTs)), 20, 'r', 'filled')
                        plot(1:length(px), px/max(abs(px)))
                        plot(1:length(vel), vel/max(vel))
                        plot(1:length(acc), acc/max(acc))
                        %                             plot(linspace(1, 10e4, length(jerk2)), jerk2/max(jerk2))
                        legend('jumpTime', 'landTime', 'px', 'vel', 'acc')
                        scatter(t, zscore(accz(t)), 20, 'b', 'filled')
                        hold off
                        title('x')
                        
                        % plot sensor input
                        subplot(2,2,3)
                        plot(-dists_x(k, jump_ts{k}(i)-120:jump_ts{k}(i)+480))
                        hold on
                        plot(-dists_y(k, jump_ts{k}(i)-120:jump_ts{k}(i)+480))
                        scatter(121, -dists_x(k, jump_ts{k}(i)), 10, 'r')
                        scatter(121, -dists_y(k, jump_ts{k}(i)), 10, 'g')
                        hold off
                        
                        subplot(2,2,4)
                        
                        fprintf([num2str(size(jumpTs)) '\n'])
                        plot(pz)
                        hold on
                        plot(jumpTs*ones(1,2), [min(pz) max(pz)], 'r') ;
                        hold off
                        shg
                        pause; %s{k}(i) = 'y';
                        prompt = input('Throw out: [Y]','s');
                        close
                        if isempty(prompt)
                            % convert to pos indices
                            maxImpLoc = jumpTs+jump_ts{k}(i)-120;
                            maxLandLoc = landTs + jump_ts{k}(i)-120;
                            %                                 maxLandLoc = round((sensorLandTs-20e3)/167)+jump_ts{k}(i);
                            
                            buffer = 600;
                            % use these values when the rat does not pause
                            % that long at each spot
                            % xPos = movmean(pos(maxImpLoc-600:maxImpLoc+300,8), 5);
                            % yPos = movmean(pos(maxImpLoc-600:maxImpLoc+300,10), 5);
                            
                            % use these values when the rat pauses a long
                            % time
                            xPos = movmean(posX(maxImpLoc-buffer:maxImpLoc+300), [0 5]);
                            yPos = movmean(posY(maxImpLoc-buffer:maxImpLoc+300), [0 5]);
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
                                if min(posX(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer)) > (min(endx)-30) & max(posX(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer)) < (max(endx)+30)
                                    
                                    % Then check and store data
                                    % If z-range is below 500
                                    if max(posZ(maxImpLoc-400:maxImpLoc+400)) - min(posZ(maxImpLoc-400:maxImpLoc+400)) < 500
                                        behavior.events.trials{count}.x =pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,8);
                                        behavior.events.trials{count}.y =pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,10);
                                        behavior.events.trials{count}.z =pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,9);
                                        behavior.events.trials{count}.mapping = nan;
                                        behavior.events.trials{count}.timestamps = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,1);
                                        behavior.events.trials{count}.errorPerMarker = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,11);
                                        
                                        
                                        %             [r1 r2 r3] = quat2angle([pos(maxImpLoc-400:maxImpLoc+400,4);pos(maxImpLoc-400:maxImpLoc+400,5)...
                                        %                 pos(maxImpLoc-400:maxImpLoc+400,6);pos(maxImpLoc-400:maxImpLoc+400,7)])
                                        behavior.events.trials{count}.orientation.x = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,4);
                                        behavior.events.trials{count}.orientation.y = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,5);
                                        behavior.events.trials{count}.orientation.z = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,6);
                                        behavior.events.trials{count}.orientation.w = pos(start_point+maxImpLoc-buffer:end_point+maxImpLoc-buffer,7);
                                        
                                        % for each jump, check which direction the rat ran
                                        %                                     if sum(diff(behavior.events.trials{count}.x)) > 0 %&& sum(diff(behavior.events.trials{count}.x)) > 0
                                        %                                         mod = 1;
                                        %                                     else
                                        %                                         mod = 2;
                                        %                                     end
                                        
                                        behavior.events.trialConditions(count) = k*2-(2-mod);% to make unique k %k * mod;
                                        behavior.events.trialIntervals(count,:) = pos([start_point+maxImpLoc-buffer end_point+maxImpLoc-buffer],1);
                                        behavior.events.jumpTime(count,:) = [pos(maxImpLoc,1) pos(maxLandLoc,1)];
                                        behavior.events.jumpLoc(count, :) = [pos(maxImpLoc, 8) pos(maxLandLoc, 8) pos(maxImpLoc, 10) pos(maxLandLoc, 10) pos(maxImpLoc, 9) pos(maxLandLoc, 9)];
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


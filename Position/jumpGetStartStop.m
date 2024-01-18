function [trialStartEnd, other] = jumpGetStartStop(behavior)
% for each trial, get the turnaround points.

trialStart = nan(length(behavior.events.trials), 1);
trialEnd = nan(length(behavior.events.trials), 1);

other = [];
figure; %hold on;
plot(behavior.position.x, behavior.position.y)
% for ii = 1:length(behavior.events.trials)
%     plot(behavior.events.trials{ii}.x, 'k')
% end
title('Click 1 = low number, click 2 = high number')

[x, y] = ginput(2);
% y(1) = -400;
% y(2) = 840;
close

for ii = 1:length(behavior.events.trials)
    direction = diff(behavior.events.trials{ii}.x([1 end])) > 0; % want to be 1 if from L->R, -1 if from R->L
    if direction < 1
        direction = -1;
    end
    
     % smooth trial a bit
    xPos = movmean(behavior.events.trials{ii}.x, 5);
    yPos = movmean(behavior.events.trials{ii}.y, 5);
    
    b = 100;
    if direction > 0
        [pks1, start_point] = min(abs(xPos(b:end-b)-x(1))+abs(yPos(b:end-b)-y(1)));
        [pks2, end_point] = min(abs(xPos(b:end-b)-x(2))+ abs(yPos(b:end-b)-y(2)));
        start_point = start_point + b;
        end_point = end_point + b;
        if pks1 > 30 || pks2 > 30
            [pks1, start_point] = min(abs(xPos-x(1))+abs(yPos-y(1)));
        	[pks2, end_point] = min(abs(xPos-x(2))+ abs(yPos-y(2)));
        end
        
    elseif direction < 0
        [pks1, start_point] = min(abs(xPos(b:end-b)-x(2))+ abs(yPos(b:end-b)-y(2)));
        [pks2, end_point] = min(abs(xPos(b:end-b)-x(1))+abs(yPos(b:end-b)-y(1)));
        start_point = start_point + b;
        end_point = end_point + b;
         if pks1 > 30 || pks2 > 30
            [pks1, start_point] = min(abs(xPos-x(2))+abs(yPos-y(2)));
        	[pks2, end_point] = min(abs(xPos-x(1))+ abs(yPos-y(1)));
        end
        
    end
    
    trialStart(ii) = behavior.events.trials{ii}.timestamps(start_point);

%     end_point = length(xPos);
%     while diff(xPos(end_point-1:end_point))*direction < 0
%         end_point = end_point-1;
%     end

    trialEnd(ii) = behavior.events.trials{ii}.timestamps(end_point);
    
    other = [other; start_point end_point];
    
    
end

trialStartEnd = [trialStart trialEnd];
% behavior.events.trialStartEndPoints = [trialStart trialEnd];
    
%% test

figure; hold on; 
for ii = 1:length(other)
    plot(movmean(behavior.events.trials{ii}.x, 5)); 
end
for ii = 1:length(other)
    scatter(other(ii, 1), behavior.events.trials{ii}.x(other(ii, 1))); 
end

for ii =1:length(other)
    scatter(other(ii,2), behavior.events.trials{ii}.x(other(ii,2))); 
end

figure; hold on; 
for ii = 1:length(other)
    plot(movmean(behavior.events.trials{ii}.y, 5)); 
end
for ii = 1:length(other)
    scatter(other(ii, 1), behavior.events.trials{ii}.y(other(ii, 1))); 
end

for ii =1:length(other)
    scatter(other(ii,2), behavior.events.trials{ii}.y(other(ii,2))); 
end



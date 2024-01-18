% convert wheel .behav file to behavior file
% tried on 3_15 DT7 so far

dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

% load([basename '.pos']);
% 
% behavior = pos2behav(pos,'behavType','jump');


load behav.mat

behavior.position.x = pos(:,8);
behavior.position.y = pos(:,10);
behavior.position.z = pos(:,9);
behavior.timestamps = pos(:,1);

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



plotCondTrials(behavior, 'dim', "xy")
plotCondTrials(behavior, 'dim', "x")
plotCondTrials(behavior, 'dim', "xz")

plotAllTrials(behavior, 'dim', "xzt")

%% somehow linearize the trials
% linearize for each trial separately
% also each condition was done twice, 5-20 min apart

% ncond = length(unique(behavior.events.trialConditions));
ncond = 4;
startx = zeros(ncond, 1); starty = zeros(ncond, 1);
endx = zeros(ncond, 1); endy = zeros(ncond, 1);

% get a starting point
for ii = 1:ncond
    ids = find(behavior.events.trialConditions == ii);
    
    
    for jj = ids
        startx(ii) = startx(ii) + behavior.events.trials{jj}.x(1);
        starty(ii) = starty(ii) + behavior.events.trials{jj}.y(1);
        
        endx(ii) = endx(ii) + behavior.events.trials{jj}.x(end);
        endy(ii) = endy(ii) + behavior.events.trials{jj}.y(end);
    end
    startx(ii) = startx(ii)/length(ids);
    starty(ii) = starty(ii)/length(ids);
    endx(ii) = endx(ii)/length(ids);
    endy(ii) = endy(ii)/length(ids);
    
end
   
behavior.events.startEndPos = [startx starty endx endy];

%% Fit a circle to the trials

% get a long vector of all x and y positions
x = []; y = [];
trials = find(behavior.events.trialConditions < 5);
for ii = trials%1:length(behavior.events.trials)
    x = [x; behavior.events.trials{ii}.x];
    y = [y; behavior.events.trials{ii}.y];
end

% now fit circle
data = [x y];

% using equations for the circle,2xc_1 + 2yc_2 + (r^2 -c_1^2-c_2^2) = x^2 +
% y^2
% and c_3 = r^1-c_1^2-c_2^2
% want to minimize the distance to this circle.
% need to determin parameters phi_1, phi_2, ..., z_1, z_2, r

A = [2*x 2*y ones(length(x), 1)];
b = x.^2 + y.^2;

% c(1) = x offset, c(2) = y offset, c(3) = r^2-c1^2 -c2^2
c = A\b;
r = sqrt(c(3)+c(1).^2+c(2).^2);

% test it out
xx = x-c(1);
yy = y-c(2);
theta = 0:pi/50:2*pi;
figure
scatter(xx, yy, 10, 'filled')
hold on
scatter(0, 0, 20, 'r', 'filled')
plot([0, 0], r*[0, 1], 'r', 'lineWidth', 1)
plot(r*cos(theta), r*sin(theta), 'k', 'lineWidth', 1)
axis square
xlabel('x')
ylabel('y')

% sweet it works

% okay so for each trial then, 
% zero the position
% get the angle of each point theta = (arctan(x/y))
% and then points along circle are [rcos(theta), rsin(theta)]

% then linearize
% get set 0 to be the start point of the trial

% then the distance along the arc is the distance along the trial.
%phew done okay.


%% Then plot position vs phase

% make bins
dx = 10;
posBins = 0:dx:2000;

dPhi = 0.1;
phaseBins = -pi:dPhi:pi;



%% test figures

figure
hold on
for ii = 1:length(behavior.events.trials)
    
    plot(behavior.events.trials{ii}.x, behavior.events.trials{ii}.y, 'k')
    axis square
end

for ii = 1:length(behavior.events.trials)
    cond = behavior.events.trialConditions(ii);
    scatter(behavior.events.startEndPos(cond, 1), behavior.events.startEndPos(cond, 2), 80, 'r', '.');
    scatter(behavior.events.startEndPos(cond, 3), behavior.events.startEndPos(cond, 4), 80, 'g', '.');
    
end

% figure
% hold on
% for ii = 1:length(behavior.events.trials)
%     
%     scatter(behavior.events.trialConditions(ii), behavior.events.trials{ii}.x(1), 20, 'k', '.')
% %     plot(behavior.events.trials{ii}.x, behavior.events.trials{ii}.y, 'k')
%     
% end
    
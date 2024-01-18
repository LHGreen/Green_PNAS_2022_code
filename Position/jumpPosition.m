% get linearized jump position
function [position, jumpLoc, landLoc] = jumpPosition(behavior)

jumpLoc = nan(13, 2);
for cond = 1:6
    trials = find(behavior.events.trialConditions == cond);
    jumpLoc(cond, 1) = nanmedian(behavior.events.jumpLoc(trials, 1));
    jumpLoc(cond, 2) = nanmedian(behavior.events.jumpLoc(trials, 3));
    jumpLoc(cond+6, 1) = nanmedian(behavior.events.jumpLoc(trials, 2));
    jumpLoc(cond+6, 2) = nanmedian(behavior.events.jumpLoc(trials, 4));
end
if length(behavior.events.startEndPos) < 3
    jumpLoc(end, 1) = behavior.events.startEndPos(1, 1); 
    jumpLoc(end, 2) = behavior.events.startEndPos(1, 2); 
else
    jumpLoc(end, 1) = behavior.events.startEndPos(1, 3);
    jumpLoc(end, 2) = behavior.events.startEndPos(1, 4);
end


% shift so that min is at 0
[position, linPoints] = linearize(behavior, 'points', jumpLoc);
landLoc = linPoints(7:12, :);
jumpLoc = linPoints(1:6, :);
position = position-linPoints(end);
landLoc = landLoc-linPoints(end);
jumpLoc = jumpLoc-linPoints(end);
jumpLoc(7:8) = nan;
landLoc(7:8) = nan;
end
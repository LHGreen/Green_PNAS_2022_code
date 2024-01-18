function [position, linPoints] = linearize(varargin)

% linearize position by projecting onto a line formed by starts and end of
% track

%%% Inputs %%%
%
% behavior:  buzcode-format behavior
%
%
% points:    additional xy points to linearize, in nx2 vector, x then y
%
%
%%% Outputs %%%
%
% position:  linearized position for entire session
%
% linPoints:    specified points to convert (if not specified, returns NaN)


p = inputParser();
addRequired(p, 'behavior', @isstruct);
addParameter(p, 'points', nan(1, 2), @isnumeric)

parse(p, varargin{:})
behavior = p.Results.behavior;
points = p.Results.points;

x = behavior.events.startEndPos(:,1);
y = behavior.events.startEndPos(:,2);

xPos = behavior.position.x-x(1);
yPos = behavior.position.y-y(1);

line = [diff(x) diff(y)]';

position = [xPos yPos]*(line/norm(line));

if sum(sum(isnan(points))) < 1
    points(:,1) = points(:,1)-x(1);
    points(:,2) = points(:,2)-y(1);
    linPoints = points*(line/norm(line));
else
    ii = isnan(points);
    points(:,1) = points(:,1)-x(1);
    points(:,2) = points(:,2)-y(1);
    linPoints = points*(line/norm(line));
end


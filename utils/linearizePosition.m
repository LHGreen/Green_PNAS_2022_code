function [position, linPoints] = linearizePosition(varargin)

% linearize position by projecting onto a line formed by starts and end of
% track

%%% Inputs %%%
%
% behavior:  buzcode-format behavior
%
% zeroPoints: [x y] coordinates of where you want the position 0 to
% correspond to

% points:    additional xy points to linearize, in nx2 vector, x then y
%
%
%%% Outputs %%%
%
% position:  linearized position for entire session
%
% linPoints:    specified points to convert (if not specified, returns NaN)
%
%
% Example:
% 
% position = linearizePosition(behavior)
% 


p = inputParser();
addRequired(p, 'behavior', @isstruct);
addParameter(p, 'zeroPoints', zeros(1, 2), @isnumeric);
addParameter(p, 'points', nan(1, 2), @isnumeric)

parse(p, varargin{:})
behavior = p.Results.behavior;
startPoints = p.Results.zeroPoints;
points = p.Results.points;

x = startPoints(1);
y = startPoints(2);

xPos = behavior.position.x-x;
yPos = behavior.position.y-y;

% do linear regression to get coefficients of line
line = [ones(length(xPos), 1) xPos]\yPos;

position = [xPos yPos]*(line/norm(line));

if sum(sum(isnan(points))) < 1
    points(:,1) = points(:,1)-x;
    points(:,2) = points(:,2)-y;
    linPoints = points*(line/norm(line));
else
    ii = isnan(points);
    points(:,1) = points(:,1)-x;
    points(:,2) = points(:,2)-y;
    linPoints = points*(line/norm(line));
end


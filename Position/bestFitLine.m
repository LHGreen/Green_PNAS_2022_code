function [angleVal] = bestFitLine(behavior, varargin)

% Get line of best fit for track, to convert x or y position to linear
% track position
% 
% INPUTS
% 
% behavior (required)
% 
% Name-value paired inputs:
% 'trials'        - Trials to include in analysis
% 
% OUTPUT
% 
% angleVal        - angle of line given by linear regression

p = inputParser;
% addRequired(p, 'behavior', @isstruct)
addParameter(p, 'trials', 1:length(behavior.events.trials), @isnumeric)

parse(p, varargin{:})
% behavior = p.Results.behavior;
trials = p.Results.trials;



% Plot all trials
trialLength = length(behavior.events.trials{1}.x);
ntrials = length(trials);
xpoints = [];
ypoints = [];

% exclude edges of trials
N = 20;

yMin = min(behavior.position.y);
xMin = min(behavior.position.x);

% yMin = 0; xMin = 0;

yLine = linspace(0, max(behavior.position.y)-yMin, 200);
xLine = linspace(0, max(behavior.position.x)-xMin, 200);

% figure
% hold on
for ii = trials
%     scatter(behavior.events.trials{ii}.x(N:end-N)-xMin, behavior.events.trials{ii}.y(N:end-N)-yMin, 15, 'k', 'filled')
    xpoints = [xpoints; behavior.events.trials{ii}.x(N:end-N)-xMin];
    ypoints = [ypoints; behavior.events.trials{ii}.y(N:end-N)-yMin];
end
b = [ones(length(xpoints), 1) xpoints]\ypoints;
% plot(xLine, b(2).*xLine + b(1), 'r', 'lineWidth', 1)

% close
angleVal = abs(atan(b(2)));


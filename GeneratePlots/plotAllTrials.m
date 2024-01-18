function [] = plotAllTrials(varargin)

p=inputParser();
addRequired(p,'behavior', @isstruct);
addParameter(p,'dim',"x",@isstring);
parse(p,varargin{:})

behavior = p.Results.behavior;
dim = p.Results.dim;

nT = length(behavior.events.trials);

figure
hold on
if dim == "x"
    for ii = 1:nT
        plot(behavior.events.trials{ii}.x)
    end
    xlabel('bin')
    ylabel('x')
end

if dim == "y"
    for ii = 1:nT
        plot(behavior.events.trials{ii}.y)
    end
    xlabel('bin')
    ylabel('y')
end

if dim == "z"
    for ii = 1:nT
        plot(behavior.events.trials{ii}.x)
    end
    xlabel('bin')
    ylabel('z')
end

if dim == "xy"
    for ii = 1:nT
        plot(behavior.events.trials{ii}.x, behavior.events.trials{ii}.y)
    end
    xlabel('x')
    ylabel('y')
end

if dim == "xz"
    for ii = 1:nT
        plot(behavior.events.trials{ii}.x, behavior.events.trials{ii}.z)
    end
    xlabel('x')
    ylabel('z')
end

if dim == "xzt"
    for ii = 1:nT
        plot3(behavior.events.trials{ii}.x, behavior.events.trials{ii}.z, behavior.events.trials{ii}.timestamps)
    end
    xlabel('x')
    ylabel('z')
    zlabel('time')
end




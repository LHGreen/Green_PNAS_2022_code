function [] = plotCondTrials(varargin)

p=inputParser();
addRequired(p,'behavior', @isstruct);
addParameter(p,'dim',"x",@isstring);
addParameter(p, 'trialType', "linear", @isstring);
parse(p,varargin{:})

behavior = p.Results.behavior;
dim = p.Results.dim;
trialType = p.Results.trialType;


ncond = unique(behavior.events.trialConditions);

figure
hold on
if dim == "x"
    for cond = ncond
        subplot(5, 2, cond)
        hold on
        trials = find(behavior.events.trialConditions == cond);
        for ii = trials
            plot(behavior.events.trials{ii}.x)
        end
    end
    xlabel('bin')
    ylabel('x')
end

if dim == "y"
    for cond = ncond
        subplot(5, 2, cond)
        hold on
        trials = find(behavior.events.trialConditions == cond);
        for ii = trials
            plot(behavior.events.trials{ii}.y)
        end
        xlabel('bin')
        ylabel('y')
    end
    
end

if dim == "z"
    for cond = ncond
        subplot(5, 2, cond)
        hold on
        trials = find(behavior.events.trialConditions == cond);
        for ii = trials
            plot(behavior.events.trials{ii}.x)
        end
        xlabel('bin')
        ylabel('z')
    end
    
end

if dim == "xy"
    for cond = ncond
        subplot(5, 2, cond)
        hold on
        trials = find(behavior.events.trialConditions == cond);
        for ii = trials
            scatter(behavior.events.trials{ii}.x, behavior.events.trials{ii}.y, 5, 'k', 'filled')
            scatter(behavior.events.trials{ii}.x(1), behavior.events.trials{ii}.y(1), 10, 'r', 'filled')
            scatter(behavior.events.trials{ii}.x(end), behavior.events.trials{ii}.y(end), 10, 'g', 'filled')
            if trialType == "jump"
                scatter(behavior.events.jumpLoc(ii, 1), behavior.events.jumpLoc(ii, 3), 10, 'c', 'filled')
                scatter(behavior.events.jumpLoc(ii, 2), behavior.events.jumpLoc(ii, 4), 10, 'm', 'filled')
            end
        end
        xlabel('x')
        ylabel('y')
    end
    
end

if dim == "xz"
    for cond = ncond
        subplot(5, 2, cond)
        hold on
        trials = find(behavior.events.trialConditions == cond);
        for ii = trials
            scatter(behavior.events.trials{ii}.x, behavior.events.trials{ii}.z, 5, 'k', 'filled')
            scatter(behavior.events.trials{ii}.x(1), behavior.events.trials{ii}.z(1), 10, 'r', 'filled')
            scatter(behavior.events.trials{ii}.x(end), behavior.events.trials{ii}.z(end), 10, 'g', 'filled')
            if trialType == "jump"
                scatter(behavior.events.jumpLoc(ii, 1), behavior.events.jumpLoc(ii, 5), 10, 'c', 'filled')
                scatter(behavior.events.jumpLoc(ii, 2), behavior.events.jumpLoc(ii, 6), 10, 'm', 'filled')
            end
        end
        xlabel('x')
        ylabel('z')
        title(['ntrials = ' num2str(length(trials))])
    end
end

if dim == "xzt"
    for cond = ncond
        subplot(5, 2, cond)
        hold on
        trials = find(behavior.events.trialConditions == cond);
        for ii = trials
            plot3(behavior.events.trials{ii}.x, behavior.events.trials{ii}.z, behavior.events.trials{ii}.timestamps)
        end
        xlabel('x')
        ylabel('z')
        zlabel('time')
    end
    
end

%% try some Kalman filter stuff

clear; close all

dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

% load behavior
load([basename '.behavior.mat']);

% But actually just do it for the z-position first
% constant acceleration model
g = 9.18*10e3/120; % 9.18 m/s
T = 1/behavior.samplingRate;

u = 0; %-g;

%%
x0 = [pos(1,3) diff(pos(1:2, 3))]';% diff(diff(pos(1:3, 3))) diff(diff(diff(pos(1:4, 3))))]';

A = [1 T ; 0 1]; % 
B = [ T 1]';

C = [1 0];% [1 0 0; 0 0 0; 0 0 0];


% A = [1 T (T^2)/2 (2*T^3)/3; 0 1 T (T^2)/2; 0 0 T 0; 0 0 0 1]; % 
% B = [(2*T^3)/3 (T^2/2) T 1]';
% 
% C = [1 0 0 0];% [1 0 0; 0 0 0; 0 0 0];

% test = [5 8 10 13 18];
% figure
% plot(yv)

% for tt = 1:length(test)
    
Q = 20; %0.1*eye(3);

R = 1; %eye(3);


sys = ss(A, B, C, 0, -1);

[Kf, L, P, M] = kalman(sys, Q, R);

% Implement the filter?
t = length(pos(:,1));

P = B*Q*B';
ye = zeros(t, 1);
ycov = zeros(t, 1);
yv = pos(1:t,3);
x = x0;

%
for ii = 1:t
    % measurement update
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

% end

%%
% Kf = lqe(A, Q, C, Q, R);
% want xyz position, velocity, acceleration
% (so a 1x9 vector?)

% u = 0;
% A = [0 1 0 0 0 0 0 0 0; % add velocity
%      0 0 1 0 0 0 0 0 0; % add acceleration
%      0 0 0 0 0 0 0 0 0;
%      0 0 0 0 1 0 0 0 0;
%      0 0 0 0 0 1 0 0 0;
%      0 0 0 0 0 0 0 0 0;
%      0 0 0 0 0 0 0 1 0;
%      0 0 0 0 0 0 0 0 1;
%      0 0 0 0 0 0 0 0 0];
%  
% B = zeros(9, 1);
% 
% % for A
% C = ones(9, 9);
% D = ones(9, 1);
% 
% Vd = eye(9); % position noise, who knows.
% 
% Vs = 0.5*eye(9); % sensor variance, usually
% 
% sys = ss(A, B, C, D);


Kf = lqe(A, Vd, C, Vd, Vn);



% get linearized position and jump position
jumpLoc = nan(13, 2);
for cond = 1:6
    trials = find(behavior.events.trialConditions == cond);
    jumpLoc(cond, 1) = median(behavior.events.jumpLoc(trials, 1));
    jumpLoc(cond, 2) = median(behavior.events.jumpLoc(trials, 3));
    jumpLoc(cond+6, 1) = median(behavior.events.jumpLoc(trials, 2));
    jumpLoc(cond+6, 2) = median(behavior.events.jumpLoc(trials, 4));
end
jumpLoc(end, 1) = behavior.events.startEndPos(1, 3);
jumpLoc(end, 2) = behavior.events.startEndPos(1, 4);

% shift so that min is at 0
[position, linPoints] = linearize(behavior, 'points', jumpLoc);
landLoc = linPoints(7:12, :);
jumpLoc = linPoints(1:6, :);
position = position-linPoints(end);
landLoc = landLoc-linPoints(end);
jumpLoc = jumpLoc-linPoints(end);



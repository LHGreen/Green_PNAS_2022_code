function [ye, ycov, R, Q] = adaptKalman(A, B, u, C, Q, R, alpha, yv, adapt)

% Kalman filter
% Constant Velocity Model
% Using Adaptive Adjustment of Noise Covariance in Kalman Filter for
% Dynamic Estimation - Akhlagi, Zhou, & Huang, Arxiv

% Inputs
% A : Transition matrix
% B : Input transition matrix
% u : forcing
% C : Output transition matrix
% Q : process noise covariance
% R : measurement noise covariance
% alpha: learning rate of adaptive adjustment process
% yv : input positions
% adapt: 1 if want to use adaptive process
% 
% Outputs
% ye: output position
% ycov: output covariance
%


P = eye(size(A));
x = [yv(1) 0]';
t = length(yv);

ye = zeros(t, 1);
ycov = zeros(t, 1);


% implement algorithm
for ii = 1:t
    % measurement update
    if sum(isnan(x)) > 0
        x = zeros(2, 1);
    end
    
    % Step 1: Prediction update
    x = A*x + B*u;
    P = A*P*A' + Q;
    
    % Step 2: Correction
    d = yv(ii)-C*x; % measurement innovation
    S = C*P*C' + R; % Innovation Covariance
    K = P*C'/S; % Kalman gain matrix
    x = x + K*d;
    P_m = P; % save for next step
    P = (eye(2) -K*C)*P; %*(eye(2)-Mn*C)' + Mn*R*Mn';
    
    ye(ii) = C*x;
    
    ycov(ii) = C*P*C';
    
    % estimate R and Q
    if adapt==1
        resid = yv(ii)-C*x;
        S_hat = resid*resid';
        R = alpha*R + (1-alpha)*(S_hat + C*P_m*C');
        
        Q = alpha*Q + (1-alpha)*(K*(d*d')*K');
    
    end
    
end



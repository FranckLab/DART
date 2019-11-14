function [x] = remove_rigid_motion(x0, x1, sizeI)
% Remove rigid body motion

% Remove points near the border of the image
idx = x0(:,1)>50 & x0(:,2)>50 & x0(:,3)>30 ...
    & x0(:,1)<sizeI(1)-50 & x0(:,2)<sizeI(2)-50 & x0(:,3)<sizeI(3)-30;

y0 = x0(idx,:);
y1 = x1(idx,:);

% Setup the optimization problem to remove the rigid drift and rotatin
% about x-y axis. 
% x1 = x0 + drift + R*x1 + deformation
% R = [m, -n, 0; n, m, 0; 0, 0, 1]; %Considering rotation only about xy
drift = optimvar('drift', 1, 3, 'LowerBound', 0, 'UpperBound', 0);
m = optimvar('m', 'LowerBound', -1, 'UpperBound', 1);
n = optimvar('n', 'LowerBound', -1, 'UpperBound', 1);

% Setup the optimization problem
prob = optimproblem;
Rx = [m*y0(:,1)'-n*y0(:,2)'; n*y0(:,1)' + m*y0(:,2)'; y0(:,3)'];
prob.Objective = sum(sum((y1 - drift - Rx).^2));

% Solve 
[sol,fval] = solve(prob);

x = y0;
end


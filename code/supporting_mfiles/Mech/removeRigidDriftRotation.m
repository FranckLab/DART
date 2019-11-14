function u = removeRigidDriftRotation(u, x)

x = cell2mat(x);

%%%%% Remove rigid translation
u = u-median(u);

%%%%% Remove xy rigid rotation;

% Prepare 
center_point = mean(x); % Approximate center about which rotation occurs
r = x(:, 1:2) - center_point(1:2); % Radius vector from center for each x
r_mag = sqrt(sum(r.^2, 2)); % Radius for each x
r_unit = r./r_mag; % Unit radius vector for each x
rotatn_mat = [0, -1; 1, 0]; % Rotatn matrix to get tangential vector from radius vector
t_unit = (rotatn_mat*(r_unit'))'; % Tangential unit vector for each x
t_mag = sum(u(:, 1:2).*t_unit, 2); % Tangential displacement magnitude

% ""Rotation angle""
theta = median(t_mag./r_mag);
t_correction = -theta*t_unit.*r_mag;

% Correct for rotation
u(:,1:2) = u(:,1:2) + t_correction;

%%%%% Remove rigid stage motion with z_height

% % % Visualize
% % % % figure; scatter(u(:,1), x(:,3))
% % % % figure; scatter(u(:,2), x(:,3))
zbucket = linspace(min(x(:,3)), max(x(:,3)), 8);
u_stage = zeros(length(zbucket)-1, 3); 
for i = 1:length(zbucket)-1
    idx = x(:,3) > zbucket(i) & x(:,3) <= zbucket(i+1);
    z(i,1) = median(x(idx, 3));
    u_stage(i,:) = median(u(idx, :));
end

% Fit spline to each displacement
u_drift_correction = zeros(size(x));
for i = 1:3
    curve = fit(z, u_stage(:,i), 'smoothingspline', 'SmoothingParam', 1);
    u_drift_correction(:,i) = -curve(x(:,3));
end
u = u + u_drift_correction;

end
% Not so simple double pendulum solved with Matlab's RK4
% ode45 solver The program simulates the double
% pendulum and the dispacements and angular velocities. 
% The bob is a point with all the mass concentrated to there
% but it is depicted as a circle, which does not have a 
% complicated behavior such as deformation

clear;
clc;

% Define the length, mass
% of the suspended point
% in SI units
L = [1 1];
M = [1 1];

% Maximum simulation time (in sec)
T_max = 20;

% Gravity acceleration g (m/s^2)
g = 9.81;

% Initial angle
theta = [3*pi/4 3*pi/4];
% Angular velocity
omega = [0.0 0.0];

% Center of the ball and radius
xc = [0 0];
r = 0.1;

% Number of points to form
% the outline of the ball
N = 2^6 + 1;

% Get circle and segment size
[x, y] = circle(xc(1),xc(2),r,N);

% Limits of the box for ploting
xmin = - sum(L) - r;
xmax = sum(L) + r;
ymin = - r;
ymax = 2 * sum(L) + r;

% Define the system of ODEs
% Inputs: time and v(angular velocity 1st pendulum,
% angular velocity 2nd pendulum,angle 1st pendulum,
% angle 2nd pendulum)
f = @(t,v) [v(3); v(4); (-g*(2*M(1) + M(2))*sin(v(1)) - M(2)*g*sin(v(1)-2*v(2))-2*sin(v(1)-v(2))*M(2)*(v(4)^2 * L(2) + v(3)^2 * L(1) * cos(v(1)-v(2))))/(L(1)*(2*M(1) + M(2) - M(2)*cos(2*v(1) - 2*v(2))));....
    (2*sin(v(1)-v(2))*(v(3)^2 * L(1) * (M(1) + M(2)) + g*(M(1) + M(2))*cos(v(1)) + v(4)^2 * L(2) * M(2) * cos(v(1)-v(2))))/(L(2)*(2*M(1) + M(2) - M(2)*cos(2*v(1) - 2*v(2))))];

% Options for the MATLAB ode45 solver
options=odeset('RelTol',1e-5,'AbsTol',1e-5);

% Solve using RK4 (ode45)
[t,vs] = ode45(f,[0 T_max],[theta,omega],options);

% Coordinate transformation
xcm_1 = [L(1)*sin(vs(:,1)), -L(1)*cos(vs(:,1))+sum(L)];
xcm_2 = [L(1)*sin(vs(:,1))+ L(2)*sin(vs(:,2)), -L(1)*cos(vs(:,1))-L(2)*cos(vs(:,2))+sum(L)];

% Angular velocities
omega_1 = vs(:,1);
omega_2 = vs(:,2);

for i = 1:size(t,1)
    time = t(i);
    % Draw everything in a single line
    xs=[0;xcm_1(i,1);x+xcm_1(i,1);xcm_1(i,1);xcm_2(i,1);x+xcm_2(i,1)];
    ys=[sum(L);xcm_1(i,2);y+xcm_1(i,2);xcm_1(i,2);xcm_2(i,2);y+xcm_2(i,2)];
    subplot(2,2,1)
    plot(xs,ys);xlabel('x (m)');ylabel('y (m)')
    % Put bounding box
    axis([xmin, xmax, ymin, ymax]);
    % Make axis equal and remove numbering
    axis square;
    
    % Getframe is required to make plotting smoother
    M=getframe;
end

% Plot separate for performance reasons
% Plot theta angles
subplot(2,2,2)
hold
plot(t,vs(:,3),'r');
plot(t,vs(:,4),'g');xlabel('t (s)');ylabel('\theta (rad)');
legend('\theta_1','\theta_2');
axis([0 T_max min([vs(:,3);vs(:,4)]) max([vs(:,3);vs(:,4)])]);
hold off    

% Plot x,y of the first center of mass
subplot(2,2,3)
hold
plot(t,xcm_1(:,1),'r');
plot(t,xcm_1(:,2),'g');xlabel('t (s)');ylabel('1st Position (m)');
legend('x(t)','y(t)');
axis([0 T_max min(xmin,ymin) max(xmax,ymax)]);
hold off

% Plot x,y of the second center of mass
subplot(2,2,4)
hold
plot(t,xcm_2(:,1),'r');
plot(t,xcm_2(:,2),'g');xlabel('t (s)');ylabel('2nd Position (m)');
legend('x(t)','y(t)');
axis([0 T_max min(xmin,ymin) max(xmax,ymax)]);
hold off


% Function that creates a circle of radius r
% and center [xc,yc]
% Input:
%   xc     : x-coordinate of center
%   yc     : y-coordinate of center
%   r      : radius of circle
%   N      : Number of points in the circle
% Output:
%   x      : x-coordinates of all points on the circle
%   y      : y-coordinate of all point on the circle
function [x, y] = circle(xc,yc,r,N)

% If radius or number of points are below prescribed
% return point
if (r <= 0) || (N <= 1)
    x = xc;
    y = yc;
    disp('Warning: r<=0 or N<=1 so only the center point is returned');
end

% Create all angles avoiding the last overlap
% of the last point
th = linspace(0,2*pi,N)';
th(end) = [];

% Create coordinates of points based on Polar
% coordinate system
x = r * cos(th) + xc;
y = r * sin(th) + yc;

end

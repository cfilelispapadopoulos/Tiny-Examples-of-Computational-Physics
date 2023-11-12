% Not so simple not so aerodynamic pendulum with a bob, with its mass
% concentrated at its center. The bob is subject to drag due to its movement
% inside the surrounding air.
% Equation: d^2 \theta /dt^2 = -g/l sin\theta - l /(2*m) \rho C_d A
% (d \theta / dt)^2 \hat{omega}
% The direction of drag force is opposite to to movement

clear;
clc;

% Define the length, mass
% of the suspended point
% in SI units
L = 1;
M = 1;
% Drag coefficient
C_d = 0.47;
% Mass density of air (kg/m^3)
rho = 1.293;

% Maximum simulation time (in sec)
T_max = 4;
t_step = 0.01;

% Define the acceleration
f = inline('-g/L*sin(x)-sign(om)*L/(2*M)*rho*C_d*A*om^2','g','L','M','rho','C_d','A','x','om');

% Gravity acceleration g (m/s^2)
g = 9.81;

% Initial angle
theta = 0.8;
% Angular velocity
omega = 0.0;

% Center of the ball and radius
xc = [0 0];
r = 0.1;
% Cross sectional area of sphere
A = pi ^ r^2;

% Center of mass
xcm = [L*sin(theta) L-L*cos(theta)];

% Number of points to form
% the outline of the ball
N = 2^6 + 1;

% Get circle and segment size
[x, y] = circle(xc(1),xc(2),r,N);

% Group xy
xy = [x-xc(1) y-xc(2)];

% Velocity of center of mass
vcm = [0 0];

% Acceleration of center of 
acm = [0 0];

% Limits of the box for ploting
xmin = - L - r;
xmax = L + r;
ymin = - r;
ymax = 2 * L + r;

% Time advancement loop
time = 0;
while time <= T_max
    % Kick (1/2)
    omega = omega + f(g,L,M,rho,C_d,A,theta,omega) * t_step / 2.0; 
    % Drift
    theta = theta + omega * t_step;
    % Convert angle to cartesian coordinates
    xcm = [L*sin(theta) L-L*cos(theta)];
    % Kick (1/2)
    omega = omega + f(g,L,M,rho,C_d,A,theta,omega) * t_step / 2.0;
    
    % Plot taking into account the position of the center of mass
    xs = [0; xcm(1)];
    xs = [xs;[xy(:,1);xy(1,1)]+xcm(1)];
    ys = [L; xcm(2)];
    ys = [ys;[xy(:,2);xy(1,2)]+xcm(2)];
    
    subplot(2,2,1)
    plot(xs,ys);
    % Put bounding box
    axis([xmin, xmax, ymin, ymax]);
    % Make axis equal and remove numbering
    axis square;

    subplot(2,2,2)
    hold
    plot(time,theta,'.r');xlabel('t (s)');ylabel('\theta (rad)');
    axis([0 T_max -pi pi]);
    hold off    
    
    subplot(2,2,3)
    hold
    plot(time,xcm(1),'.r');xlabel('t (s)');ylabel('x (m)');
    axis([0 T_max xmin xmax]);
    hold off
    
    subplot(2,2,4)
    hold
    plot(time,xcm(2),'.r');xlabel('t (s)');ylabel('y (m)');
    axis([0 T_max ymin ymax]);
    hold off    
    
    % Advance time variable
    time = time + t_step;
    % Getframe is required to make plotting smoother
    MM=getframe;

end

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

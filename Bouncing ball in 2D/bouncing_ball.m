% Simulation of a circle (ball)
% the resistance to deformation 
% is simulated with
% a series of springs

clear;
clc;

% Define prescribed values and constants
% Center of box
xb = 0.5;
yb = 0.5;

% Center of circle
xc = 0.5;
yc = 0.4;
xcm = [xc yc];

% Radius
r = 0.1;

% Define number of points in circle N
N = 2^6 + 1;

% Mass per point of the ball
mass = 1 / (N - 1);

% Gravity acceleration g (m/s^2) of circle
g = 9.81;

% Spring constants for rebound of the ball
% surface
spr_coeff = 60;
% Critical damping coefficient
spr_damp = sqrt(4*mass*spr_coeff);

% Define maximum time of simulation in seconds
% and timestep
T_max = 1.5;
t_step = 0.001;
t_step_min = t_step;

% Get circle and segment size
[x, y] = circle(xc,yc,r,N);

% Group xy
xy = [x-xc y-yc];
xy0 = xy;

% Velocity of points
vxy = zeros(size(xy));
% Velocity of center of mass
vcm = [0 2];

% Acceleration of points
axy = zeros(size(vxy));
% Acceleration of center of 
acm = [0 0];

% Define box of the simulation as 10 times the radius
scale_box = 10;

% Maximum width and height of the box
w_max = scale_box * r;
h_max = scale_box * r;

% Limits of the box for ploting
xmin = xb - w_max / 2;
xmax = xb + w_max / 2;
ymin = yb - h_max / 2;
ymax = yb + h_max / 2;

% Time advancement loop
time = 0;
while time <= T_max
    % Kick (1/2)
    vxy = vxy + axy * t_step / 2.0;
    vcm = vcm + acm * t_step / 2.0; 
    % Drift
    xy = xy + vxy * t_step;
    xcm = xcm + vcm * t_step;
    % Check positions and correct velocities and positions
    [xy, vxy] = checkPos(xy, vxy, xcm, xb, yb, w_max, h_max);
    % update accelerations
    [axy, acm] = computeAcc(xy, xy0, vxy, spr_coeff, spr_damp, g, r, mass);
    % Kick (1/2)
    vxy = vxy + axy * t_step / 2.0;
    vcm = vcm + acm * t_step / 2.0;
    
    % Plot taking into account the position of the center of mass
    plot([xy(:,1)+xcm(1);xy(1,1)+xcm(1);xcm(1)],[xy(:,2)+xcm(2);xy(1,2)+xcm(2);xcm(2)]);
    
    % Put bounding box
    axis([xmin, xmax, ymin, ymax]);
    % Make axis equal and remove numbering
    axis equal;
    axis off;
    
    % Advance time variable
    time = time + t_step;
    % Getframe is required to make plotting smoother
    M=getframe;

end

% Function that computes accelerations
% Inputs:
%   xy_          : Positions
%   xy0_         : Spring rest positions
%   vxy_         : Velocities
%   spr_coeff    : Spring coefficient
%   damp_coeff   : Damping coefficient
%   g_           : Gravity acceleration
%   r_           : Radius of the ball
%   mass_        : Total mass of circle
% Outputs:
%   axy_         : Accelerations of all points
%   acm_         : Acceleration of center of mass
function [axy_, acm_] = computeAcc(xy_, xy0_, vxy_, spr_coeff_, spr_damp_, g_, r_, mass_)
    % Matrix to retain accelerations per dimension
    axy_ = zeros(size(xy_));
    acm_ = [0 0];
    
    % Add gravity acceleration towards -y direction
    acm_(2) = - g_;
    
    % Length of springs
    l_xy = sqrt(sum(xy_.^2,2));
    % Spring displacement from equilibrium
    dxy = (xy_-xy0_);
    ind = (l_xy < (1-1e-3)*r_) | (l_xy > (1+1e-3)*r_);
    
    % Acceleration due to spring contraction or elongation
    gacc = -spr_coeff_/mass_*dxy-spr_damp_/mass_*vxy_;
    
    % Assign acceleration to points
    axy_(ind,:)=gacc(ind,:);
    
    % Add reaction forces to center of mass
    gacccm=sum(gacc(ind,:));
    acm_=acm_-gacccm/(mass_*size(axy_,1));
end

% Function that checks for collisions and flips positions
% and velocieties
% Inputs:
%   xy_    : Positions of points
%   vxy_   : Velocities of points
%   xcm_   : Position of CM
%   xb_    : X-coordinate of center
%   yb_    : Y-coordinate of center
%   w_max_ : Width of box
%   h_max_ : Heigth of box
% Outputs:
%   xy_    : New positions of points
%   vxy_   : New velocities of points
function [xy_, vxy_] = checkPos(xy_, vxy_, xcm_, xb_, yb_, w_max_, h_max_)
    % Coordinates of the up, down, left, right most points
    y_up = yb_ + h_max_ / 2;
    y_down = yb_ - h_max_ / 2;
    x_left = xb_ - w_max_ / 2;
    x_right = xb_ + w_max_ / 2;
    
    % Check up / down
    % Up
    ind = find(xy_(:,2) + xcm_(2) > y_up);
    xy_(ind,2) = 2 * y_up - xy_(ind,2) - 2*xcm_(2);
    vxy_(ind,2) = 0;
        
    % Down
    ind = find(xy_(:,2) + xcm_(2) < y_down);
    xy_(ind,2) = 2* y_down - xy_(ind,2) - 2*xcm_(2);
    vxy_(ind,2) = 0;
    
    % Check left / right
    % Left
    ind = find(xy_(:,1) + xcm_(1) < x_left);
    xy_(ind,1) = 2 * x_left - xy_(ind,1) - 2*xcm_(1);
    vxy_(ind,1) = 0;
    
    % Right
    ind = find(xy_(:,1) + xcm_(1) > x_right);
    xy_(ind,1) = 2 * x_right - xy_(ind,1) - 2*xcm_(1);
    vxy_(ind,1) = 0;
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

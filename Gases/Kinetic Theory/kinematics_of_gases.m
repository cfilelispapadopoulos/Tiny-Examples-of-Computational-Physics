% Oxygen filled piston in vacuum compressed
% with a spring. The molecules are considered as
% groups with velocities obtained from Maxwell - Boltzmann
% Distribution. This groups (particles) are not point
% They are modelled as squares with side length dictated by the density 
% inside the cylinder \rho = n * molar_mass / (height * piston_position)
% Initially everything is in equilibrium
% Inner pressure is atmospheric pressure
% The spring is at rest
%
% Authors note: Not very accurate simulation but good for qualitative
% showcase

clear;
clc;
% All units are in SI
% Width of piston
w = 0.3;
% Height of piston
h = 0.1;
% Thickness of piston
t = 0.01;
% Mass of piston
Mt = 0.2;
% Initial position of the piston
p = 0.2;
pin = p;
% Safeguard against choosing an initial position that towards the very end
% of the available length of the piston
p = min(p, 0.75 * w);

% Spring coeff
spr_coeff = 3;

% Maximum simulation time
T_max = 10;
% Time step
t_step = 0.005;

% Temperature of Gas in K
T = 500;
% Boltzmann's constant (m^2 kg s^-2 K^-1)
kB = 1.380649e-23;
% Ideal Gas Constant (J/(K mol))
R = 8.314;
% Avogadro's number
NA = 6.0221408e+23;
% Number of moles
n = 20;
% Number of particles (in simulation)
N = 1000;
% Mass of each particle (in simulation)
% 1 mol of O_2 is 32g (molar mass)
mmol = 32e-3;
M = n * mmol / N;


% Create initial positions and velocities
[pp,vp] = create_initial_pos_vel(N, h, p, kB, T, mmol, NA);
% Accelerations of particles
ap = zeros(size(pp));

% Velocity of wall
vw = 0;
% Acceleration of wall
aw = 0;

% Create a figure
figure
% Loop over
% Initialize time
time = 0;
while time < T_max
    % Clear previous plot
    clf
    % Compute adaptive time step
    t_step=min(0.001,0.01/max(max((abs(vp)))));
    % Hold the plots in order to simultaneously present the geometry
    % and particles
    hold
    % Cleate geometry
    create_geometry(w, h, t, p);

    % Leapfrog KDK (particles)
    % Kick (1/2)
    vp = vp + ap * t_step / 2.0; 
    % Drift
    pp = pp + vp * t_step;
    % Kick (1/2)
    vp = vp + ap * t_step / 2.0;

    % Leapfrog KDK (wall)
    % Kick (1/2)
    vw = vw + aw * t_step / 2.0; 
    % Drift
    p = p + vw * t_step;

    % Compute new density
    rho = mmol * n / (h * p);
    % Kick (1/2)
    vw = vw + aw * t_step / 2.0;

    % Check position with respect to the box and 
    % compute acceleration of the lid
    [pp, vp, aw] = checkPosWF(pp, vp, h, p, M, Mt, rho);
    % Compute the acceleration due to the spring
    aw = aw + cylinderAcc(spr_coeff, pin, p, Mt);
    % Plot particles
    plot(pp(:,1), pp(:,2),'.g');

    % Update time
    time = time + t_step;
    % Added for visualization reasons
    MM=getframe;
    % Remove hold to move to next frame
    hold off
end

% Function that creates the geometry of piston and spring
% Inputs:
%   w_  : Width of piston
%   h_  : Height of piston
%   t_  : Thickness of piston
%   p_  : Position of piston
function create_geometry(w_, h_, t_, p_)
    % Cylinder
    plot([0 w_], [0 0], 'k', 'LineWidth', 3);
    plot([0 w_], [h_ h_], 'k', 'LineWidth', 3);
    plot([0 0], [0 h_], 'k', 'LineWidth', 3);
    plot([w_], [h_/2],'ok');
    % Piston
    plot([p_ p_], [0 h_], 'r', 'LineWidth', 3); 
    plot([p_ + t_ p_ + t_], [0 h_], 'r', 'LineWidth', 3);
    axis([-0.1, w_+0.1, -0.1, h_+0.1]);
    % Plot spring using the following lib
    % https://www.mathworks.com/matlabcentral/fileexchange/28724-2d-spring-coordinates-for-plotting
    [x, y] = springcoord([p_ + t_, h_/2], [w_, h_/2], 5, w_, 0.05*(w_-(p_+t_)));
    % Plot spring
    plot(x,y);
end

% Function that creates the initial positions using uniform random numbers
% and velocities following Maxwell - Boltzmann distribution
% Inputs:
%   N_     : Number of particles (groups)
%   h_     : Height of the piston
%   p_     : Position of the piston with respect to reference point
%   kB_    : Boltzmann constant
%   T_     : Temperature (K)
%   mmol_  : Molar mass of O_2
%   NA_    : Avogadro's number
% Outputs:
%   pos_   : Position of particles
%   vel_   : Velocities of particles
function [pos_, vel_] = create_initial_pos_vel(N_, h_, p_, kB_, T_, mmol_, NA_)
    % Create position of particles in the cylinder
    pos_ = rand(N_,2).*[p_, h_];
    v = sqrt((kB_ * T_) / (mmol_ / NA_));
    vel_ = v * randn(N_, 2);
end

% Function that checks position of particles and compute force on collision
% on the piston
% Inputs:
%   pos_:   Position of all particles
%   cel_:   Velocities of all particles
%   h_  :   Height of piston
%   p_  :   Position of piston relative to refernce (0,0)
%   M_  :   Mass of each square particle
%   Mt_ :   Mass of the piston
%   rho_:   Average density inside the piston
% Outputs:
%   pos_:   New positions after collision
%   vel_:   New velocities after collision
%   aw_ :   Acceleration due to collision
function [pos_, vel_, aw_] = checkPosWF(pos_, vel_, h_, p_, M_, Mt_, rho_)
    N = size(pos_,1);
    % Check and reflect (elastic collisions)
    % Left
    ind = (pos_(:,1) < 0);
    pos_(ind,1) = - pos_(ind,1);
    vel_(ind,1) = - vel_(ind,1);
    
    % Right
    ind = (pos_(:,1) > p_);
    pos_(ind,1) = 2 * p_ - pos_(ind,1);
    % Sum before reflection
    % Account for the difference in time
    % to collition between the particles in the 
    % group (box of side d)
    d = sqrt(M_ / (rho_));
    % Force spread across the square (F dx dy \approx F d d)
    aw_ = M_ / Mt_ * sum((vel_(ind,1).^2) / (p_)) * d^2;
    % Reflect velocity across x axis
    vel_(ind,1) = - vel_(ind,1);
    
    % Down
    ind = (pos_(:,2) < 0);
    pos_(ind,2) = - pos_(ind,2);
    vel_(ind,2) = - vel_(ind,2);
    
    % Up
    ind = (pos_(:,2) > h_);
    pos_(ind,2) = 2 * h_- pos_(ind,2);
    vel_(ind,2) = - vel_(ind,2);    
end

% Function to compute acceleration of the piston due to the spring
% Inputs:
%   spr_coeff : Spring coeffient
%   pin_      : Initial position of the piston
%   p_        : Current position of the piston
%   Mt_       : Mass of the piston
% Outputs:
%   aw_       : Acceleration of the piston
function [aw_] = cylinderAcc(spr_coeff_, pin_, p_, Mt_)
    aw_ = -spr_coeff_ * (p_ - pin_)/ Mt_;
end